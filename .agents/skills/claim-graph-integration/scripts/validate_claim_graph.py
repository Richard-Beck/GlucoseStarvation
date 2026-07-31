#!/usr/bin/env python3
"""Validate the clean claim graph schema."""

import argparse
import json
import sys
from pathlib import Path


SCHEMA_VERSION = "claim-graph/v4"
REQUIRED_TOP_LEVEL = ("metadata", "claims", "evidence", "relationships")
RELATIONS = {"support", "undermine"}
STRENGTHS = {"strong", "moderate", "weak"}
FORBIDDEN_CLAIM_FIELDS = {
    "qualified_by_claims",
    "qualifies",
    "supported_by_claims",
    "supports",
    "undermines",
    "evidence",
}
def require_list(value, path, errors):
    if isinstance(value, list):
        return value
    errors.append(f"{path} must be a list")
    return []


def unique_nodes(nodes, section, errors):
    result = {}
    for index, node in enumerate(nodes):
        if not isinstance(node, dict):
            errors.append(f"{section}[{index}] must be an object")
            continue
        node_id = str(node.get("id", "")).strip()
        if not node_id:
            errors.append(f"{section}[{index}] is missing id")
        elif node_id in result:
            errors.append(f"{section} has duplicate id: {node_id}")
        else:
            result[node_id] = node
    return result


def validate_authoritative_contract(metadata, claims, errors):
    contract = metadata.get("authoritative_claim_contract")
    if not isinstance(contract, dict):
        errors.append("metadata.authoritative_claim_contract must be an object")
        return
    fixed_claims = require_list(
        contract.get("claims"),
        "metadata.authoritative_claim_contract.claims",
        errors,
    )
    fixed_ids = set()
    for index, fixed in enumerate(fixed_claims):
        path = f"metadata.authoritative_claim_contract.claims[{index}]"
        if not isinstance(fixed, dict):
            errors.append(f"{path} must be an object")
            continue
        fixed_id = str(fixed.get("id", "")).strip()
        fixed_text = str(fixed.get("text", "")).strip()
        if not fixed_id or not fixed_text:
            errors.append(f"{path} requires id and text")
            continue
        if fixed_id in fixed_ids:
            errors.append(f"authoritative contract has duplicate id: {fixed_id}")
            continue
        fixed_ids.add(fixed_id)
        claim = claims.get(fixed_id)
        if claim is None:
            errors.append(f"authoritative claim is missing: {fixed_id}")
            continue
        if claim.get("user_fixed") is not True:
            errors.append(f"authoritative claim {fixed_id} must have user_fixed=true")
        for field, expected in fixed.items():
            if claim.get(field) != expected:
                errors.append(
                    f"authoritative claim {fixed_id}.{field} differs from contract: "
                    f"{claim.get(field)!r} != {expected!r}"
                )
    for claim_id, claim in claims.items():
        if not isinstance(claim.get("user_fixed"), bool):
            errors.append(f"claim {claim_id}.user_fixed must be boolean")
        elif claim.get("user_fixed") and claim_id not in fixed_ids:
            errors.append(
                f"claim {claim_id} has user_fixed=true but is absent from the "
                "authoritative contract"
            )


def validate_graph(graph):
    errors = []
    if not isinstance(graph, dict):
        return ["graph must be an object"]
    for key in REQUIRED_TOP_LEVEL:
        if key not in graph:
            errors.append(f"missing top-level key: {key}")

    metadata = graph.get("metadata")
    if not isinstance(metadata, dict):
        errors.append("metadata must be an object")
        metadata = {}
    if metadata.get("schema_version") != SCHEMA_VERSION:
        errors.append(
            "metadata.schema_version must be "
            f"{SCHEMA_VERSION!r}, got {metadata.get('schema_version')!r}"
        )
    relation_values = set(
        require_list(metadata.get("relation_values"), "metadata.relation_values", errors)
    )
    strength_values = set(
        require_list(metadata.get("strength_values"), "metadata.strength_values", errors)
    )
    if relation_values != RELATIONS:
        errors.append(
            "metadata.relation_values must contain only support and undermine"
        )
    if strength_values != STRENGTHS:
        errors.append(
            "metadata.strength_values must contain strong, moderate, and weak"
        )

    claim_nodes = require_list(graph.get("claims"), "claims", errors)
    evidence_nodes = require_list(graph.get("evidence"), "evidence", errors)
    relationship_nodes = require_list(
        graph.get("relationships"),
        "relationships",
        errors,
    )
    claims = unique_nodes(claim_nodes, "claims", errors)
    evidence = unique_nodes(evidence_nodes, "evidence", errors)
    relationships = unique_nodes(relationship_nodes, "relationships", errors)

    for claim_id, claim in claims.items():
        if not str(claim.get("text", "")).strip():
            errors.append(f"claim {claim_id} is missing text")
        forbidden = sorted(FORBIDDEN_CLAIM_FIELDS.intersection(claim))
        if forbidden:
            errors.append(
                f"claim {claim_id} uses obsolete embedded relationship fields: "
                f"{', '.join(forbidden)}"
            )
    for evidence_id, node in evidence.items():
        if not str(node.get("figure_id", "")).strip():
            errors.append(f"evidence {evidence_id} is missing figure_id")
        panels = require_list(node.get("panels"), f"evidence {evidence_id}.panels", errors)
        if not panels or any(not str(panel).strip() for panel in panels):
            errors.append(f"evidence {evidence_id}.panels must be nonempty strings")
        if not str(node.get("observation", "")).strip():
            errors.append(f"evidence {evidence_id} is missing observation")
        if not str(node.get("source", "")).strip():
            errors.append(f"evidence {evidence_id} is missing source")

    used_evidence = set()
    for relationship_id, relationship in relationships.items():
        source_type = relationship.get("source_type")
        source_id = str(relationship.get("source_id", "")).strip()
        target_id = str(relationship.get("target_claim_id", "")).strip()
        relation = relationship.get("relation")
        strength = relationship.get("strength")
        reason = str(relationship.get("reason", "")).strip()

        if source_type not in {"evidence", "claim"}:
            errors.append(
                f"relationship {relationship_id}.source_type must be evidence or claim"
            )
        elif source_type == "evidence":
            if source_id not in evidence:
                errors.append(
                    f"relationship {relationship_id} references unknown evidence "
                    f"source: {source_id}"
                )
            else:
                used_evidence.add(source_id)
        elif source_id not in claims:
            errors.append(
                f"relationship {relationship_id} references unknown claim source: "
                f"{source_id}"
            )

        if target_id not in claims:
            errors.append(
                f"relationship {relationship_id} references unknown target claim: "
                f"{target_id}"
            )
        if source_type == "claim" and source_id and source_id == target_id:
            errors.append(f"relationship {relationship_id} is a claim self-edge")
        if relation not in RELATIONS:
            errors.append(
                f"relationship {relationship_id}.relation must be support or "
                "undermine"
            )
        if strength not in STRENGTHS:
            errors.append(
                f"relationship {relationship_id}.strength must be strong, "
                "moderate, or weak"
            )
        if not reason:
            errors.append(f"relationship {relationship_id} is missing reason")

    unused = sorted(set(evidence).difference(used_evidence))
    if unused:
        errors.append(
            "evidence nodes without relationships: " + ", ".join(unused)
        )

    validate_authoritative_contract(metadata, claims, errors)
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", type=Path)
    args = parser.parse_args()
    try:
        graph = json.loads(args.graph.read_text(encoding="utf-8"))
    except Exception as exc:
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
        f"{len(graph.get('relationships', []))} relationships"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
