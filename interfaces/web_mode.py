from __future__ import annotations
import os
import shlex
import traceback
import tempfile
from pathlib import Path
from flask import Flask, request, jsonify, send_from_directory



def run_web(
    *,
    args,
    collect_linkers,
    collect_nodes,
    assemble,
    run_xtb_for_out_dir,
    run_obabel_for_out_dir,
    normalize_shapes
) -> None:

    ROOT_DIR = Path(__file__).resolve().parent.parent
    STATIC_DIR = ROOT_DIR / "static"
    
    app = Flask(__name__, static_folder=str(STATIC_DIR))

    ### Store dynamic paths in Flask's config
    app.config['PATHS'] = {
        'linkers': None,
        'nodes': None,
        'out': None,
        'assemblies': None
    }
    
    shapes_list = normalize_shapes(args.shape)



    ### Routes

    @app.route("/")
    def index():
        return app.send_static_file("index.html")

    @app.route("/api/defaults", methods=["GET"])
    def get_defaults():
        """Sends the CLI arguments to the frontend to pre-fill the path inputs."""
        return jsonify({
            "nodes": str(Path(args.nodes).resolve()) if args.nodes else "",
            "linkers": str(Path(args.linkers).resolve()) if args.linkers else "",
            "out": str(Path(args.out).resolve()) if args.out else ""
        })

    @app.route("/api/config", methods=["POST"])
    def set_config():
        data = request.json
        if not data:
            return jsonify({"error": "No data provided"}), 400

        l_dir = Path(data.get("linkers", "")).expanduser().resolve()
        n_dir = Path(data.get("nodes", "")).expanduser().resolve()
        o_dir = Path(data.get("out", "")).expanduser().resolve()

        if not l_dir.is_dir():
            return jsonify({"error": f"Linkers directory not found: {l_dir}"}), 400
        if not n_dir.is_dir():
            return jsonify({"error": f"Nodes directory not found: {n_dir}"}), 400

        try:
            o_dir.mkdir(parents=True, exist_ok=True)
            a_dir = o_dir / "assemblies"
            a_dir.mkdir(parents=True, exist_ok=True)
            
            app.config['PATHS']['linkers'] = l_dir
            app.config['PATHS']['nodes'] = n_dir
            app.config['PATHS']['out'] = o_dir
            app.config['PATHS']['assemblies'] = a_dir

            ### Fetch linkers/nodes
            linkers = [f.name for f in l_dir.glob("*") if f.is_file() and f.suffix in {".xyz", ".mol", ".json"}]
            nodes = [f.name for f in n_dir.glob("*.json") if f.is_file()]
            
            ### Scan for previously built cages
            existing_cages = []
            if a_dir.exists():
                for d in sorted(a_dir.iterdir(), key=lambda x: x.stat().st_mtime, reverse=True):
                    if d.is_dir():
                        xtb_file = d / "xtb_opt.xyz"
                        ff_file = d / "ff_opt.mol"
                        base_file = d / f"{d.name}.mol"
                        
                        if xtb_file.exists() and xtb_file.stat().st_size > 0:
                            existing_cages.append({
                                "name": d.name,
                                "url": f"/files/assemblies/{d.name}/{xtb_file.name}",
                                "tag": "Optimized (XTB)"
                            })
                        elif ff_file.exists() and ff_file.stat().st_size > 0:
                            existing_cages.append({
                                "name": d.name,
                                "url": f"/files/assemblies/{d.name}/{ff_file.name}",
                                "tag": "Optimized (OBABEL)"
                            })
                        elif base_file.exists() and base_file.stat().st_size > 0:
                            existing_cages.append({
                                "name": d.name,
                                "url": f"/files/assemblies/{d.name}/{base_file.name}",
                                "tag": "Assembled"
                            })

            return jsonify({
                "message": "Directories loaded successfully.",
                "linkers": sorted(linkers),
                "nodes": sorted(nodes),
                "shapes": shapes_list,
                "history": existing_cages
            })
        except Exception as e:
            traceback.print_exc()
            return jsonify({"error": str(e)}), 500

    @app.route("/api/build", methods=["POST"])
    def build_cage():
        """Triggers the assembly of a specific node and linker."""
        if not app.config['PATHS']['linkers']:
            return jsonify({"error": "Please set directories first."}), 400

        data = request.json
        if not data:
            return jsonify({"error": "No JSON payload provided"}), 400

        node_name = data.get("node")
        linker_name = data.get("linker")
        shape = data.get("shape")

        if not all([node_name, linker_name, shape]):
            return jsonify({"error": "Missing node, linker, or shape"}), 400

        L_DIR = app.config['PATHS']['linkers']
        N_DIR = app.config['PATHS']['nodes']
        O_DIR = app.config['PATHS']['out']
        A_DIR = app.config['PATHS']['assemblies']

        try:
            ### Converts .xyz/.mol to .json if necessary
            linker_jsons = collect_linkers(L_DIR, O_DIR)
            node_jsons = collect_nodes(N_DIR, O_DIR)

            ### Find the exact JSON files needed for this build
            target_node = next((n for n in node_jsons if n.stem == Path(node_name).stem), None)
            target_linker = next((l for l in linker_jsons if l.stem == Path(linker_name).stem), None)

            if not target_node or not target_linker:
                return jsonify({"error": "Could not locate preprocessed JSON for selected node/linker"}), 404

            ### Assemble
            assemble(target_node, target_linker, A_DIR, shape)
            
            ### Move the .mol file into its own subfolder
            for f in A_DIR.glob("[!.]*"):
                if f.is_file() and f.suffix == ".mol":
                    dest = A_DIR / f.stem
                    dest.mkdir(parents=True, exist_ok=True)
                    f.replace(dest / f.name)
            
            ### Identify the newly created directory to send back to the frontend
            prefix = f"{target_node.stem}__{target_linker.stem}__"
            built_dirs = [d for d in A_DIR.glob(f"{prefix}*") if d.is_dir()]
            
            if not built_dirs:
                return jsonify({"error": "Assembly completed but output directory not found."}), 500
                
            ### Grab the most recently modified directory matching the prefix
            built_dir = max(built_dirs, key=lambda d: d.stat().st_mtime)
            built_mol = built_dir / f"{built_dir.name}.mol"
            
            if not built_mol.exists():
                return jsonify({"error": f"Directory {built_dir.name} created, but no .mol found."}), 500

            return jsonify({
                "message": "Assembly successful",
                "output_dir": built_dir.name,
                "file_url": f"/files/assemblies/{built_dir.name}/{built_mol.name}"
            })

        except Exception as e:
            traceback.print_exc()
            return jsonify({"error": str(e)}), 500

    @app.route("/api/optimize", methods=["POST"])
    def optimize_cage():
        """Triggers xTB or OpenBabel optimization on a generated cage."""
        if not app.config['PATHS']['assemblies']:
            return jsonify({"error": "Please set directories first."}), 400

        data = request.json
        if not data:
            return jsonify({"error": "No JSON payload provided"}), 400

        target_dir_name = data.get("target_dir")
        method = data.get("method", "xtb")
        flags_str = data.get("flags", "")

        if not target_dir_name:
            return jsonify({"error": "Missing target_dir"}), 400

        A_DIR = app.config['PATHS']['assemblies']
        target_dir = A_DIR / target_dir_name
        if not target_dir.is_dir():
            return jsonify({"error": "Target directory not found"}), 404

        extra_flags = shlex.split(flags_str) if flags_str else []

        try:
            with tempfile.TemporaryDirectory() as tmp:
                tmp_parent = Path(tmp)
                tmp_target = tmp_parent / target_dir.name
                tmp_target.symlink_to(target_dir, target_is_directory=True)
                
                if method == "xtb":
                    run_xtb_for_out_dir(tmp_parent, extra_flags=extra_flags)
                    result_mol = target_dir / "xtb_opt.xyz"
                    failed_mol = target_dir / "xtb_failed_last.xyz"
                elif method == "obabel":
                    run_obabel_for_out_dir(tmp_parent, extra_flags=extra_flags)
                    result_mol = target_dir / "ff_opt.mol"
                    failed_mol = target_dir / "ff_failed_last.mol"
                else:
                    return jsonify({"error": f"Unknown optimization method: {method}"}), 400

            ### Check if successful
            if result_mol.exists() and result_mol.stat().st_size > 0:
                return jsonify({
                    "message": f"Optimization ({method}) successful",
                    "file_url": f"/files/assemblies/{target_dir.name}/{result_mol.name}"
                })
            elif failed_mol.exists():
                return jsonify({
                    "error": f"Optimization failed. Output written to {failed_mol.name}.",
                    "file_url": f"/files/assemblies/{target_dir.name}/{failed_mol.name}"
                }), 500
            else:
                return jsonify({"error": "Optimization failed to produce any output."}), 500

        except Exception as e:
            traceback.print_exc()
            return jsonify({"error": str(e)}), 500

    @app.route("/files/<path:filepath>")
    def serve_files(filepath):
        """Serves generated JSON, XYZ, and MOL files to the frontend viewer."""
        O_DIR = app.config['PATHS']['out']
        if not O_DIR:
            return "Output directory not configured", 404
        return send_from_directory(O_DIR, filepath)



    ### Server Startup

    port = int(os.environ.get("PORT", 5001))
    print(f"\n[INFO] Starting Cageinator Web Server on http://0.0.0.0:{port}\n")
    
    ### Run the Flask app (using host="0.0.0.0" to expose it outside the Docker container)
    app.run(host="0.0.0.0", port=port, debug=False)
