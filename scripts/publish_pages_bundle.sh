#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUNDLE_DIR=${1:-"$PROJECT_ROOT/latest_manuscript_bundle"}
RELEASE_TAG=${RELEASE_TAG:-pages-site}
RELEASE_ASSET=${RELEASE_ASSET:-manuscript-pages.tar.gz}
WORKFLOW=${WORKFLOW:-publish-report.yml}

if [[ ! -d "$BUNDLE_DIR" ]]; then
  echo "Bundle directory does not exist: $BUNDLE_DIR" >&2
  exit 1
fi

if [[ ! -s "$BUNDLE_DIR/index.html" ]]; then
  echo "Bundle must contain a non-empty index.html: $BUNDLE_DIR" >&2
  exit 1
fi

for command in gh tar; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "Required command is unavailable: $command" >&2
    exit 1
  fi
done

if find "$BUNDLE_DIR" -type l -print -quit | grep -q .; then
  echo "Bundle must not contain symbolic links" >&2
  exit 1
fi

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/manuscript-pages.XXXXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT
archive="$tmp_dir/$RELEASE_ASSET"

tar --hard-dereference -czf "$archive" -C "$BUNDLE_DIR" .

if gh release view "$RELEASE_TAG" >/dev/null 2>&1; then
  gh release upload "$RELEASE_TAG" "$archive" --clobber
else
  gh release create "$RELEASE_TAG" "$archive" \
    --title "GitHub Pages manuscript bundle" \
    --notes "Mutable static-site bundle used by the Pages deployment workflow." \
    --prerelease
fi

gh workflow run "$WORKFLOW"

echo "Uploaded $BUNDLE_DIR and triggered $WORKFLOW"
