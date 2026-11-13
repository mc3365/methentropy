"""Command line interface: ``methentropy run | aggregate | calibrate | lookup``."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

from . import __version__


def _workers(value):
    if value is not None:
        return value
    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    return int(slurm) if slurm else 1


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="methentropy",
        description="Depth-matched null models for 4-CpG methylation entropy (RLM output).")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="compare every window with its depth-matched references")
    r.add_argument("-i", "--input", required=True,
                   help="RLM entropy output (23 columns, optionally .gz); legacy 5-column "
                        "files need --legacy-coverage")
    r.add_argument("-o", "--output", required=True, help="output TSV (.gz to compress)")
    r.add_argument("--nsims", type=int, default=1000, help="Monte Carlo draws per null key")
    r.add_argument("--seed", type=int, default=1, help="random seed (reproducible output)")
    r.add_argument("--workers", type=int, help="processes (default: $SLURM_CPUS_PER_TASK or 1)")
    r.add_argument("--min-coverage", type=int, default=10,
                   help="minimum reads per window (RLM default is 10)")
    r.add_argument("--max-coverage", type=int, default=1000,
                   help="skip windows with more reads (repeats/artefacts)")
    r.add_argument("--chunksize", type=int, default=1_000_000)
    r.add_argument("--with-pooled", action="store_true",
                   help="also report the v1/v2 pooled null at the window's depth")
    r.add_argument("--legacy-coverage", type=int,
                   help="assumed depth for legacy 5-column input (not depth-matched)")

    a = sub.add_parser("aggregate", help="summarise windows over BED regions")
    a.add_argument("-i", "--input", required=True, help="output of 'methentropy run'")
    a.add_argument("-r", "--regions", required=True, help="BED: chrom start end [name]")
    a.add_argument("-o", "--output", required=True)
    a.add_argument("--q", type=float, default=0.05, help="FDR threshold for window calls")

    c = sub.add_parser("calibrate", help="check the test on linkage-free synthetic data")
    c.add_argument("-i", "--input", required=True, help="RLM entropy output (23 columns)")
    c.add_argument("--windows", type=int, default=2000)
    c.add_argument("--nsims", type=int, default=1000)
    c.add_argument("--seed", type=int, default=1)
    c.add_argument("--min-coverage", type=int, default=10)
    c.add_argument("--max-coverage", type=int, default=1000)
    c.add_argument("--details", help="optional TSV with per-window calibration results")

    lk = sub.add_parser("lookup", help="write exact minimum / bimodal entropy for given depths")
    lk.add_argument("--nreads", type=int, nargs="+", required=True)
    lk.add_argument("-o", "--output", required=True)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", stream=sys.stderr)

    if args.command == "run":
        from .pipeline import run
        run(args.input, args.output, n_sims=args.nsims, seed=args.seed,
            workers=_workers(args.workers), min_coverage=args.min_coverage,
            max_coverage=args.max_coverage, chunksize=args.chunksize,
            with_pooled=args.with_pooled, legacy_coverage=args.legacy_coverage)
    elif args.command == "aggregate":
        from .aggregate import aggregate
        aggregate(args.input, args.regions, args.output, q_threshold=args.q)
    elif args.command == "calibrate":
        from .calibrate import calibrate
        summary, details = calibrate(args.input, n_windows=args.windows, n_sims=args.nsims,
                                     seed=args.seed, min_coverage=args.min_coverage,
                                     max_coverage=args.max_coverage)
        print(json.dumps(summary, indent=2))
        if args.details:
            details.to_csv(args.details, sep="\t", index=False, float_format="%.6g")
    elif args.command == "lookup":
        import pandas as pd
        from .core import N_CPGS, bimodal_entropy, minimum_entropy
        rows = [(n, k, k / (N_CPGS * n), minimum_entropy(n, k), bimodal_entropy(n, k))
                for n in args.nreads for k in range(N_CPGS * n + 1)]
        pd.DataFrame(rows, columns=["n_reads", "meth_sites", "methylation", "entropy_min",
                                    "entropy_bimodal"]).to_csv(args.output, sep="\t",
                                                               index=False, float_format="%.6g")
    return 0


def entry():
    try:
        return main()
    except (ValueError, FileNotFoundError) as err:
        print(f"methentropy: error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(entry())
