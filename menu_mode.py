
from __future__ import annotations
import shlex
import time
from itertools import product
from pathlib import Path



##########################################################################
### NOTICE: This feature is deprecated and will be removed in the future
##########################################################################



def run_menu(
    *,
    collect_linkers,
    collect_nodes,
    assemble,
    run_xtb_for_out_dir,
    SHAPES: dict[str, object]) -> None:

    xtb_opt_enabled = False
    xtb_flags: list[str] = []
    out_root: Path | None = None
    nodes_dir: Path | None = None
    linkers_dir: Path | None = None
    linkers_preprocessed = False
    chosen_shapes: list[str] | None = None

    def sep() -> None:
        print("\n" + "=" * 64 + "\n")

    print("\nWelcome to the Cageinator")
    time.sleep(1)

    while True:
        sep()
        print("******** Main Menu ********\n")
        print("1. Set output destination folder")
        print("2. Select folder with Node JSON files")
        print("3. Select folder with Linker XYZ files")
        print("4. (Optional) Pre-create framework from Linker XYZ to JSON")
        print("5. (Optional) Choose geometric shape (empty = all)")
        print("6. Toggle optimization of all built cages using xTB")
        print("7. (Optional) Set additional xTB-opt flags")
        print("8. Start assembly")
        print("9. Exit\n")
        try:
            choice = input("Enter your choice (1-9): ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\n\nExiting...\n")
            return

        if choice == "1":
            p = Path(input("\nDestination Path: ").strip()).expanduser().resolve()
            p.mkdir(parents=True, exist_ok=True)
            out_root = p
            print(f"\nOutput set to: {out_root}")

        elif choice == "2":
            p = Path(input("\nNode-Folder Path: ").strip()).expanduser().resolve()
            if not p.is_dir():
                print("\nInvalid path.")
                input("Press Enter...")
                continue
            nodes_dir = p
            print(f"\nNodes: {nodes_dir}")

        elif choice == "3":
            p = Path(input("\nLinker-Folder Path: ").strip()).expanduser().resolve()
            if not p.is_dir():
                print("\nInvalid path.")
                input("Press Enter...")
                continue
            linkers_dir = p
            print(f"\nLinkers: {linkers_dir}")

        elif choice == "4":
            if not out_root or not linkers_dir or not nodes_dir:
                print("\nSet destination (1), node folder (2), and linker folder (3) first.")
                input("Press Enter...")
                continue
            try:
                produced_linkers = collect_linkers(linkers_dir, out_root)
                produced_nodes = collect_nodes(nodes_dir, out_root)
                linkers_preprocessed = True if produced_linkers else False
                print(f"\nCreated/updated {len(produced_linkers)} linker JSON files at {out_root/'linkers_json'}")
                print(f"Created/updated {len(produced_nodes)} node JSON files at {out_root/'nodes_json'}")
            except Exception as e:
                print(f"\nError: {e}")
            input("\nPress Enter to continue...")

        elif choice == "5":
            raw = input("\nShapes (comma-separated). Leave empty for all: ").strip()
            if not raw:
                chosen_shapes = None
                print("\nShapes set to: ALL")
            else:
                toks = [t.strip() for t in raw.split(",") if t.strip()]
                invalid = [t for t in toks if t not in SHAPES]
                if invalid:
                    print(f"\nInvalid shapes: {', '.join(invalid)}")
                    input("Press Enter...")
                    continue
                chosen_shapes = toks
                print("\nShapes set to: " + ", ".join(chosen_shapes))

        elif choice == "6":
            xtb_opt_enabled = not xtb_opt_enabled
            print(f"\nxTB optimization: {'ON' if xtb_opt_enabled else 'OFF'}")

        elif choice == "7":
            raw = input(
                "\nEnter xtb_opt.sh flags (quoted, space-separated). Empty clears.\n"
                "Examples: --method gfnff --threads 8 --namespace job42\n> "
            ).strip()
            if not raw:
                xtb_flags = []
                print("\nxTB flags cleared.")
            else:
                try:
                    xtb_flags = shlex.split(raw)
                    print("\nxTB flags set to: " + " ".join(xtb_flags))
                except ValueError as e:
                    print(f"\nInvalid flags: {e}")
            input("\nPress Enter to continue...")

        elif choice == "8":
            if not out_root or not nodes_dir or not linkers_dir:
                print("\nSet destination (1), nodes (2), and linkers (3) first.")
                input("Press Enter...")
                continue
            try:
                if not linkers_preprocessed:
                    print("\nNo pre-created framework detected. Creating now...\n")
                    produced_linkers = collect_linkers(linkers_dir, out_root)
                    produced_nodes = collect_nodes(nodes_dir, out_root)
                    linkers_preprocessed = True if produced_linkers else False
                    print(f"\nCreated/updated {len(produced_linkers)} linker JSON files.")
                    print(f"\nCreated/updated {len(produced_nodes)} node JSON files.\n")

                shapes = chosen_shapes if chosen_shapes else list(SHAPES.keys())
                linker_json_dir = out_root / "linkers_json"
                nodes_json_dir = out_root / "nodes_json"
                linker_json = sorted(linker_json_dir.glob("*.json"))
                node_json = sorted(nodes_json_dir.glob("*.json"))

                if not linker_json:
                    print("\nNo linker JSONs found after preprocessing.")
                    input("Press Enter...")
                    continue
                if not node_json:
                    print("\nNo node JSONs found after preprocessing.")
                    input("Press Enter...")
                    continue

                for nj, lj in product(node_json, linker_json):
                    for shape in shapes:
                        out_dir = out_root / "assemblies"
                        assemble(nj, lj, out_dir, shape)

                        for f in out_dir.glob("[!.]*"):
                            if not f.is_file():
                                continue
                            dest = out_dir / f.stem
                            dest.mkdir(parents=True, exist_ok=True)
                            f.replace(dest / f.name)

                if xtb_opt_enabled:
                    ### Pass through any extra flags if the implementation supports it
                    try:
                        run_xtb_for_out_dir(out_dir, extra_flags=xtb_flags)     ### type: ignore[arg-type]
                    except TypeError:
                        run_xtb_for_out_dir(out_dir)                            ### fallback if function has no extra_flags

                print("\nSuccess.")

            except Exception as e:
                print(f"\nError: {e}")
            input("\nPress Enter to continue...")

        elif choice == "9":
            print("\nExiting...\n")
            break

        else:
            print("\nInvalid choice, please select 1-9.")
            input("Press Enter to continue...")
