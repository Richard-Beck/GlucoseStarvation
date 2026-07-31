# v03 A/I/D serving bundle

This package applies the ordered A/I/D bundle contract to the user-approved
v03 Results, canonical v03 Methods, and user-approved v03 figure legends.

`current_bundle/context.txt` concatenates those three sources in that order.
The valid v02 A/I/D output is preserved as `old_bundle_1`. Because the current
bundle has no newly authored `aid_source.md`, the files under `served/` are
byte-preserving heading-delimited extracts from
`old_bundle_1/aid_source.md`. No A/I/D prose was created, reviewed, or revised.

`status.json` records the source hashes, bundle order, prior-root disposition,
served source, and validation outcome. `build_bundle.py` is the deterministic
builder and refuses to overwrite an existing generated bundle.
