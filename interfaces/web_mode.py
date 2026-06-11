from __future__ import annotations
import os
import shlex
import traceback
import tempfile
from pathlib import Path
from flask import Flask, request, jsonify, send_from_directory, send_file
import re
import datetime
import shutil



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
    STATIC_DIR = Path(__file__).resolve().parent / "static"
    
    app = Flask(__name__, static_folder=str(STATIC_DIR))
    shapes_list = normalize_shapes(args.shape)

    active_sessions = {}

    def get_session_paths(session_id: str) -> dict:
        if session_id not in active_sessions:
            sys_temp = Path(tempfile.gettempdir()) / "cageinator_workspaces"
            base = sys_temp / f"session_{session_id}"
            
            l_dir = base / "linkers"
            n_dir = base / "nodes"
            o_dir = base / "output"
            a_dir = o_dir / "assemblies"
            
            l_dir.mkdir(parents=True, exist_ok=True)
            n_dir.mkdir(parents=True, exist_ok=True)
            a_dir.mkdir(parents=True, exist_ok=True)
            
            active_sessions[session_id] = {
                'base': base,
                'linkers': l_dir,
                'nodes': n_dir,
                'out': o_dir,
                'assemblies': a_dir
            }
        return active_sessions[session_id]



    ### Routes

    @app.route("/")
    def index():
        return app.send_static_file("index.html")

    @app.route("/api/sync_session", methods=["POST"])
    def sync_session():
        data = request.json
        session_id = data.get("session_id")
        if not session_id:
            return jsonify({"error": "No session ID provided"}), 400

        paths = get_session_paths(session_id)

        try:
            ### Fetch linkers/nodes
            linkers = [f.name for f in paths['linkers'].glob("*") if f.is_file() and f.suffix in {".xyz", ".mol", ".json"}]
            nodes = [f.name for f in paths['nodes'].glob("*.json") if f.is_file()]
            
            ### Scan for previously built cages
            existing_cages = []
            if paths['assemblies'].exists():
                for d in sorted(paths['assemblies'].iterdir(), key=lambda x: x.stat().st_mtime, reverse=True):
                    if d.is_dir():
                        xtb_file = d / "xtb_opt.xyz"
                        ff_file = d / "ff_opt.mol"
                        base_file = d / f"{d.name}.mol"
                        
                        if xtb_file.exists() and xtb_file.stat().st_size > 0:
                            existing_cages.append({"name": d.name, "url": f"/files/{session_id}/assemblies/{d.name}/{xtb_file.name}", "tag": "Optimized (XTB)"})
                        elif ff_file.exists() and ff_file.stat().st_size > 0:
                            existing_cages.append({"name": d.name, "url": f"/files/{session_id}/assemblies/{d.name}/{ff_file.name}", "tag": "Optimized (OBABEL)"})
                        elif base_file.exists() and base_file.stat().st_size > 0:
                            existing_cages.append({"name": d.name, "url": f"/files/{session_id}/assemblies/{d.name}/{base_file.name}", "tag": "Assembled"})

            return jsonify({
                "message": "Session loaded successfully.",
                "linkers": sorted(linkers),
                "nodes": sorted(nodes),
                "shapes": shapes_list,
                "history": existing_cages
            })
        except Exception as e:
            traceback.print_exc()
            return jsonify({"error": str(e)}), 500

    @app.route("/api/upload", methods=["POST"])
    def upload_file():
        session_id = request.form.get("session_id")
        file_type = request.form.get("type") 
        if not session_id or 'file' not in request.files:
            return jsonify({"error": "Missing file or session ID"}), 400
        
        paths = get_session_paths(session_id)
        file = request.files['file']
        
        if file.filename == '':
            return jsonify({"error": "No selected file"}), 400
            
        target_dir = paths.get(file_type)
        if not target_dir:
            return jsonify({"error": "Invalid upload type"}), 400
            
        try:
            file_path = target_dir / file.filename
            file.save(file_path)
            
            linkers = [f.name for f in paths['linkers'].glob("*") if f.is_file() and f.suffix in {".xyz", ".mol", ".json"}]
            nodes = [f.name for f in paths['nodes'].glob("*.json") if f.is_file()]
            
            return jsonify({
                "message": f"Successfully uploaded {file.filename}",
                "linkers": sorted(linkers),
                "nodes": sorted(nodes),
                "shapes": shapes_list
            })
        except Exception as e:
            traceback.print_exc()
            return jsonify({"error": str(e)}), 500

    @app.route("/api/build", methods=["POST"])
    def build_cage():
        """Triggers the assembly of a specific node and linker."""
        data = request.json
        session_id = data.get("session_id")
        node_name = data.get("node")
        linker_name = data.get("linker")
        shape = data.get("shape")

        if not all([session_id, node_name, linker_name, shape]):
            return jsonify({"error": "Missing parameters"}), 400

        paths = get_session_paths(session_id)

        try:
            ### Converts .xyz/.mol to .json if necessary
            linker_jsons = collect_linkers(paths['linkers'], paths['out'])
            node_jsons = collect_nodes(paths['nodes'], paths['out'])

            ### Find the exact JSON files needed for this build
            target_node = next((n for n in node_jsons if n.stem == Path(node_name).stem), None)
            target_linker = next((l for l in linker_jsons if l.stem == Path(linker_name).stem), None)

            if not target_node or not target_linker:
                return jsonify({"error": "Could not locate preprocessed JSON for selected node/linker"}), 404

            ### Assemble
            assemble(target_node, target_linker, paths['assemblies'], shape)
            
            newest_dir = None
            for f in paths['assemblies'].glob("[!.]*"):
                if f.is_file() and f.suffix == ".mol":
                    dest = paths['assemblies'] / f.stem
                    dest.mkdir(parents=True, exist_ok=True)
                    f.replace(dest / f.name)
                    newest_dir = dest
            
            if newest_dir:
                built_dir = newest_dir
            else:
                prefix = f"{target_node.stem}__{target_linker.stem}__"
                target_shape = shape.strip().lower()
                built_dirs = []

                search_str_1 = f"_{target_shape}__"
                search_str_2 = f"_{target_shape}"
                
                for d in paths['assemblies'].glob(f"{prefix}*"):
                    if d.is_dir():
                        dir_name = d.name.lower()
                        
                        if search_str_1 in dir_name or dir_name.endswith(search_str_2):
                            built_dirs.append(d)
                
                if not built_dirs:
                    return jsonify({"error": "Assembly completed but output directory not found."}), 500
                    
                def get_mtime(d):
                    mol_file = d / f"{d.name}.mol"
                    if mol_file.exists():
                        return mol_file.stat().st_mtime_ns
                    return d.stat().st_mtime_ns
                    
                built_dir = max(built_dirs, key=get_mtime)
                
            built_mol = built_dir / f"{built_dir.name}.mol"
            
            return jsonify({
                "message": "Assembly successful",
                "output_dir": built_dir.name,
                "file_url": f"/files/{session_id}/assemblies/{built_dir.name}/{built_mol.name}"
            })

        except Exception as e:
            traceback.print_exc()
            return jsonify({"error": str(e)}), 500

    @app.route("/api/optimize", methods=["POST"])
    def optimize_cage():
        """Triggers xTB or OpenBabel optimization on a generated cage."""
        data = request.json
        session_id = data.get("session_id")
        target_dir_name = data.get("target_dir")
        method = data.get("method", "xtb")  
        flags_str = data.get("flags", "")

        if not session_id or not target_dir_name:
            return jsonify({"error": "Missing parameters"}), 400

        paths = get_session_paths(session_id)
        target_dir = paths['assemblies'] / target_dir_name
        
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
                    "file_url": f"/files/{session_id}/assemblies/{target_dir.name}/{result_mol.name}"
                })
            elif failed_mol.exists():
                return jsonify({
                    "error": f"Optimization failed. Output written to {failed_mol.name}.",
                    "file_url": f"/files/{session_id}/assemblies/{target_dir.name}/{failed_mol.name}"
                }), 500
            else:
                return jsonify({"error": "Optimization failed to produce any output."}), 500

        except Exception as e:
            traceback.print_exc()
            return jsonify({"error": str(e)}), 500

    @app.route("/api/download", methods=["GET"])
    def download_results():
        session_id = request.args.get("session_id")
        target_dir_name = request.args.get("target_dir")
        
        if not session_id or not target_dir_name: return "Missing parameters", 400
        
        paths = get_session_paths(session_id)
        target_dir = paths['assemblies'] / target_dir_name
        
        if not target_dir.is_dir(): return "Directory not found", 404
            
        tmp_dir = Path(tempfile.gettempdir())
        zip_base_path = tmp_dir / target_dir_name
        zip_file_path = tmp_dir / f"{target_dir_name}.zip"
        shutil.make_archive(str(zip_base_path), 'zip', target_dir)
        
        return send_file(zip_file_path, as_attachment=True)

    @app.route("/api/download_workspace", methods=["GET"])
    def download_workspace():
        session_id = request.args.get("session_id")
        if not session_id: return "Missing session ID", 400
        
        paths = get_session_paths(session_id)
        tmp_dir = Path(tempfile.gettempdir())
        
        timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
        zip_filename = f"Cageinator_Workspace_{timestamp}"
        
        zip_base_path = tmp_dir / zip_filename
        zip_file_path = tmp_dir / f"{zip_filename}.zip"
        
        shutil.make_archive(str(zip_base_path), 'zip', paths['base'])
        return send_file(zip_file_path, as_attachment=True)

    @app.route("/api/clear_session", methods=["POST"])
    def clear_session():
        data = request.json
        session_id = data.get("session_id")
        if not session_id: return jsonify({"error": "Missing session ID"}), 400
        
        if session_id in active_sessions:
            paths = active_sessions[session_id]
            shutil.rmtree(paths['base'], ignore_errors=True)
            del active_sessions[session_id]
            
        get_session_paths(session_id)
            
        return jsonify({"message": "Session cleared"})

    ### Server Startup
    @app.route("/files/<session_id>/<path:filepath>")
    def serve_files(session_id, filepath):
        paths = get_session_paths(session_id)
        return send_from_directory(paths['out'], filepath)

    port = int(os.environ.get("PORT", 5001))
    print(f"\n[INFO] Starting Cageinator Web Server on http://0.0.0.0:{port}\n")

    ### Run the Flask app (using host="0.0.0.0" to expose it outside the Docker container)
    app.run(host="0.0.0.0", port=port, debug=False)
