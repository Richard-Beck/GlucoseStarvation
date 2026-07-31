# Graph Analysis Workflow

Use this optional workflow to deterministically analyze or reduce a completed
method-table provenance graph. The scripts here are for graph inspection,
ranking, and bounded work-packet construction. They do not write manuscript
Methods prose and should not be treated as automatic Methods-section generators.

Required input for node ranking:

```text
id | parent | what | why | comment
```

Required inputs for method-atom graph reduction:

```text
id | parent | what | why | comment
id | classification
```

The classification table must follow the Node classification optional workflow
and represent every provenance-table `id` exactly once.

## Workflows

### Rank provenance graph nodes

Use `scripts/rank_provenance_graph_nodes.py` when the task asks for graph
metrics, root distance, descendant counts, leaf counts, or ranking nodes by
downstream manuscript/panel-like support.

Example:

```bash
python3 .agents/skills/method-table-provenance/scripts/rank_provenance_graph_nodes.py \
  path/to/method_table.md \
  --output path/to/graph_metrics.md
```

Optional TSV output:

```bash
python3 .agents/skills/method-table-provenance/scripts/rank_provenance_graph_nodes.py \
  path/to/method_table.md \
  --format tsv \
  --output path/to/graph_metrics.tsv
```

Primary output columns include:

- `root_distance`: shortest graph distance from any parsed root row.
- `direct_children`: number of immediate downstream row IDs.
- `descendant_rows`: transitive downstream row count.
- `leaf_descendants`: downstream leaf rows.
- `manuscript_panel_like_leaf_descendants`: downstream leaf rows matching the
  script's manuscript/panel-like endpoint heuristic.

Treat manuscript/panel-like leaf detection as a heuristic. If a formal manifest
defines the panel universe for the task, compare against that manifest before
using counts as manuscript-facing evidence.

### Build method-atom graph

Use `scripts/build_method_atom_graph.py` when the task asks for deterministic
graph reduction into manageable computation-centered work packets.

Example:

```bash
python3 .agents/skills/method-table-provenance/scripts/build_method_atom_graph.py \
  path/to/method_table.md \
  path/to/node_classification.md \
  --output-dir path/to/method_atom_graph
```

The script creates one method atom per `computation_step` row. Atom core
membership is intentionally narrow: the center computation plus direct
non-computation parents and children. Longer non-computation walks are reported
as context and atom-to-atom bridge paths rather than folded into the atom core.

Main outputs:

- `method_atoms.tsv`: one row per computation-centered atom.
- `method_atom_edges.tsv`: atom-to-atom edges and bridge nodes.
- `method_atom_node_membership.tsv`: direct atom core membership.
- `method_atom_context_membership.tsv`: longer upstream/downstream context.
- `uncovered_noncomputation_nodes.tsv`: classified nodes not covered by any
  computation-centered neighborhood.
- `method_atom_proxy_candidates.tsv`: connected uncovered regions with
  manuscript/panel-like support that may need reclassification, artifact-centered
  treatment, or explicit interpretation.
- `method_atom_graph_summary.md`: compact run summary and top-ranked atoms.

## Interpretation Rules

Use graph analysis outputs as deterministic evidence for downstream reasoning,
not as manuscript prose.

Do not assume a high-ranking graph node or atom should become a Methods section.
High downstream support means the node is reused, not necessarily that it is a
scientific method.

Do not assume `computation_step` rows are inherently Methods-worthy. Many are
implementation plumbing such as reshaping, indexing, local CSV materialization,
plot-object creation, or figure-polishing wrappers.

Do not ignore `code_or_config_reference` rows. Some are disposable configs, but
some encode central scientific methods such as model definitions, likelihoods,
observation models, model-family lists, priors, or inference job specifications.
Promote these interpretively when they define methods a reader must understand.

Use method atoms to create bounded review packets or audit graph structure. For
Methods writing, build a separate interpretive layer of Methods topics/modules
that may collapse many atoms, promote important code/config nodes, and demote
display-only or audit-only branches.

Useful follow-up questions after graph analysis:

- Which high-support nodes are scientific concepts rather than file plumbing?
- Which `code_or_config_reference` nodes define model, observation, fitting, or
  simulation methods?
- Which `computation_step` atoms should be collapsed into a broader data
  assembly or analysis module?
- Which proxy candidates represent missing computation-step classifications?
- Which branches are display/layout provenance only?
