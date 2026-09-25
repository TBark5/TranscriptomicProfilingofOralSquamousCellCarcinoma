"""CLI entry points; no network or synthetic fallback on import."""
import argparse
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Paired OSCC discovery and independent validation")
    parser.add_argument("command", choices=["download", "run", "validate", "external-validation"])
    parser.add_argument("--dataset", choices=["GSE20116", "GSE184616"], default="GSE20116")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--alpha", type=float, default=.05)
    parser.add_argument("--lfc", type=float, default=1.)
    parser.add_argument("--min-count", type=int, default=10)
    parser.add_argument("--min-samples", type=int, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--libraries", nargs="+", default=["MSigDB_Hallmark_2020"],
                        choices=["MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "Reactome_2022"])
    args = parser.parse_args()
    os.environ.setdefault("MPLCONFIGDIR", str(args.root / ".tools/matplotlib"))
    try:
        if args.command == "download":
            from .acquisition import acquire
            print(acquire(args.root, args.offline, args.dataset))
        elif args.command == "validate":
            from .dashboard import load_results
            bundle = load_results(args.root, args.dataset)
            print(f"Verified complete {bundle['manifest']['dataset']} output bundle")
        else:
            from .pipeline import run
            options = vars(args).copy()
            options.pop("command")
            if args.command == "external-validation":
                options["dataset"] = "GSE184616"
            run(**options)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"Analysis error: {exc}\nNo synthetic results were substituted.\n")


if __name__ == "__main__":
    main()
