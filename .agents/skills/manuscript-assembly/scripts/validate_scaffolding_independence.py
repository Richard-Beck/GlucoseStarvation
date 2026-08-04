#!/usr/bin/env python3
"""Replay an assembly with only its declared live and runtime dependencies."""

from __future__ import print_function

import argparse
import sys
from pathlib import Path

from _assembly_contract import GateError, load_context, validate_independence


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int)
    parser.add_argument("--keep-work-root", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    context = load_context(args.project_root, args.config)
    return validate_independence(
        *context,
        timeout_override=args.timeout_seconds,
        keep_work_root=args.keep_work_root
    )


if __name__ == "__main__":
    try:
        sys.exit(main())
    except GateError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(2)
