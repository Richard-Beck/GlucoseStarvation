#!/usr/bin/env python3
"""Prepare project-local contracts for a durable manuscript assembly."""

from __future__ import print_function

import argparse
import sys
from pathlib import Path

from _assembly_contract import GateError, load_context, materialize_candidates


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--output", help="override the assembly-relative candidate manifest path"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    project_root, _, config, assembly_root = load_context(
        args.project_root, args.config
    )
    return materialize_candidates(
        config, project_root, assembly_root, output_override=args.output
    )


if __name__ == "__main__":
    try:
        sys.exit(main())
    except GateError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(2)
