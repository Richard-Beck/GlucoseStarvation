#!/usr/bin/env python3
"""Validate the manuscript claim/evidence graph prototype.

This is intentionally lightweight: it checks JSON syntax, required top-level
sections, unique node ids, and internal references from claims to claims,
evidence, methods, and assumptions.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REQUIRED_TOP_LEVEL_KEYS = ("metadata", "claims", "evidence", "methods", "assumptions")
CLAIM_REF_FIELDS = ("supported_by_claims", "qualified_by_claims", "supports", "qualifies", "undermines")
EVIDENCE_REF_FIELD = "evidence"
METHOD_REF_FIELD = "depends_on_methods"
ASSUMPTION_REF_FIELD = "depends_on_assumptions"
USER_FIXED_CONTRACT_KEY = "user_fixed_claim_contract"
USER_FIXED_LOCKED_FIELDS = ("text", "tier", "status")


def ids_by_section(graph: dict) -> dict[str, set[str]]:
    return {
        section: {node.get("id") for node in graph.get(section, [])}
        for section in ("claims", "evidence", "methods", "assumptions")
    }


def duplicate_ids(nodes: list[dict]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for node in nodes:
        node_id = node.get("id")
        if node_id in seen:
            duplicates.add(node_id)
        seen.add(node_id)
    return sorted(duplicates)


def require_list(value: object, path: str, errors: list[str]) -> list:
    if isinstance(value, list):
        return value
    errors.append(f"{path} must be a list")
    return []


def check_refs(
    refs: list[str],
    valid_ids: set[str],
    path: str,
    errors: list[str],
) -> None:
    for ref in refs:
        if ref not in valid_ids:
            errors.append(f"{path} references unknown id: {ref}")


def validate_user_fixed_contract(graph: dict, errors: list[str]) -> None:
    metadata = graph.get("metadata", {})
    contract = metadata.get(USER_FIXED_CONTRACT_KEY, {})
    contract_claims = require_list(contract.get("claims", []), f"metadata.{USER_FIXED_CONTRACT_KEY}.claims", errors)
    if not contract_claims:
        return

    claims_by_id = {claim.get("id"): claim for claim in graph.get("claims", []) if isinstance(claim, dict)}
    fixed_ids: set[str] = set()

    for index, fixed_claim in enumerate(contract_claims):
        if not isinstance(fixed_claim, dict):
            errors.append(f"metadata.{USER_FIXED_CONTRACT_KEY}.claims[{index}] must be an object")
            continue

        fixed_id = fixed_claim.get("id")
        if not fixed_id:
            errors.append(f"metadata.{USER_FIXED_CONTRACT_KEY}.claims[{index}] is missing id")
            continue
        if fixed_id in fixed_ids:
            errors.append(f"metadata.{USER_FIXED_CONTRACT_KEY}.claims has duplicate id: {fixed_id}")
        fixed_ids.add(fixed_id)

        claim = claims_by_id.get(fixed_id)
        if claim is None:
            errors.append(f"user_fixed claim is missing from claims: {fixed_id}")
            continue
        if claim.get("user_fixed") is not True:
            errors.append(f"claim {fixed_id} is listed in user_fixed contract but user_fixed is not true")

        for field in USER_FIXED_LOCKED_FIELDS:
            if field not in fixed_claim:
                errors.append(f"metadata.{USER_FIXED_CONTRACT_KEY}.claims[{fixed_id}] is missing {field}")
            elif claim.get(field) != fixed_claim.get(field):
                errors.append(
                    f"claim {fixed_id}.{field} differs from user_fixed contract: "
                    f"{claim.get(field)!r} != {fixed_claim.get(field)!r}"
                )

    for claim_id, claim in claims_by_id.items():
        if "user_fixed" in claim and not isinstance(claim.get("user_fixed"), bool):
            errors.append(f"claim {claim_id}.user_fixed must be boolean when present")
        if claim.get("user_fixed") is True and claim_id not in fixed_ids:
            errors.append(f"claim {claim_id} has user_fixed=true but is absent from metadata.{USER_FIXED_CONTRACT_KEY}")


def validate_graph(graph: dict) -> list[str]:
    errors: list[str] = []

    for key in REQUIRED_TOP_LEVEL_KEYS:
        if key not in graph:
            errors.append(f"missing top-level key: {key}")

    for section in ("claims", "evidence", "methods", "assumptions"):
        require_list(graph.get(section, []), section, errors)

    node_ids = ids_by_section(graph)
    for section in ("claims", "evidence", "methods", "assumptions"):
        for index, node in enumerate(graph.get(section, [])):
            if not isinstance(node, dict):
                errors.append(f"{section}[{index}] must be an object")
                continue
            if not node.get("id"):
                errors.append(f"{section}[{index}] is missing id")

        duplicates = duplicate_ids(graph.get(section, []))
        if duplicates:
            errors.append(f"{section} has duplicate ids: {', '.join(duplicates)}")

    metadata = graph.get("metadata", {})
    status_values = set(metadata.get("status_values", []))
    claim_tiers = set(metadata.get("claim_tiers", []))

    for index, claim in enumerate(graph.get("claims", [])):
        if not isinstance(claim, dict):
            continue
        claim_id = claim.get("id", f"claims[{index}]")
        path = f"claim {claim_id}"

        status = claim.get("status")
        if status_values and status not in status_values:
            errors.append(f"{path} has status outside metadata.status_values: {status}")

        tier = claim.get("tier")
        if claim_tiers and tier not in claim_tiers:
            errors.append(f"{path} has tier outside metadata.claim_tiers: {tier}")

        for field in CLAIM_REF_FIELDS:
            refs = require_list(claim.get(field, []), f"{path}.{field}", errors)
            check_refs(refs, node_ids["claims"], f"{path}.{field}", errors)

        evidence_refs = require_list(claim.get(EVIDENCE_REF_FIELD, []), f"{path}.{EVIDENCE_REF_FIELD}", errors)
        check_refs(evidence_refs, node_ids["evidence"], f"{path}.{EVIDENCE_REF_FIELD}", errors)

        method_refs = require_list(claim.get(METHOD_REF_FIELD, []), f"{path}.{METHOD_REF_FIELD}", errors)
        check_refs(method_refs, node_ids["methods"], f"{path}.{METHOD_REF_FIELD}", errors)

        assumption_refs = require_list(claim.get(ASSUMPTION_REF_FIELD, []), f"{path}.{ASSUMPTION_REF_FIELD}", errors)
        check_refs(assumption_refs, node_ids["assumptions"], f"{path}.{ASSUMPTION_REF_FIELD}", errors)

    validate_user_fixed_contract(graph, errors)

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "graph",
        nargs="?",
        default=Path(__file__).with_name("manuscript_claim_evidence_graph_v3_prototype.json"),
        type=Path,
        help="Path to claim/evidence graph JSON.",
    )
    args = parser.parse_args()

    try:
        graph = json.loads(args.graph.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - validator should report parse/read failures cleanly.
        print(f"ERROR: could not read JSON from {args.graph}: {exc}", file=sys.stderr)
        return 2

    errors = validate_graph(graph)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(
        "OK: "
        f"{len(graph.get('claims', []))} claims, "
        f"{len(graph.get('evidence', []))} evidence nodes, "
        f"{len(graph.get('methods', []))} methods, "
        f"{len(graph.get('assumptions', []))} assumptions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
