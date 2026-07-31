#!/usr/bin/env python3
"""Validate manuscript figure-polishing input and output contracts."""

import argparse
import csv
import json
import os
import re
import struct
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
PATH_LIKE_RE = re.compile(
    r"^[A-Za-z0-9_./ -]+\.(png|pdf|csv|tsv|txt|md|json|rds|Rds|RDS|R|Rmd)$"
)
Y_INVERSION_RE = re.compile(
    r"(1\s*-\s*[^\n;]*(?:y_npc|height_npc))|"
    r"(unit\s*\(\s*1\s*,\s*['\"]npc['\"]\s*\)\s*-\s*[^\n;]*(?:y_npc|height_npc|\by\b))"
)
RASTER_ASSEMBLY_RE = re.compile(
    r"\b(readPNG|image_read|draw_image|grid\.raster|rasterGrob)\s*\("
)
RASTER_CALL_RE = re.compile(
    r"\b(readPNG|image_read|draw_image|grid\.raster|rasterGrob)\s*\((.*?)\)",
    re.IGNORECASE | re.DOTALL,
)
PNG_LITERAL_RE = re.compile(r"['\"]([^'\"]+\.png)['\"]", re.IGNORECASE)
APPROVED_RASTER_DIRNAME = "user-approved-raster-figures"
IN_FIGURE_TITLE_RE = re.compile(
    r"\bggtitle\s*\(|"
    r"\blabs\s*\([^)]*\b(?:title|subtitle|caption)\s*=|"
    r"\bplot_annotation\s*\([^)]*\b(?:title|subtitle|caption)\s*=|"
    r"\b(draw_label|annotate)\s*\([^)]*(?:Figure\s+S?\d+|title|subtitle)",
    re.IGNORECASE | re.DOTALL,
)
BASE_PROVENANCE_COLUMNS = {
    "figure",
    "panel",
    "generator",
    "command",
    "data_inputs",
    "output_image",
    "notes",
}
LEGACY_PROVENANCE_COLUMNS = BASE_PROVENANCE_COLUMNS | {"source_image"}
REGENERATED_PROVENANCE_COLUMNS = BASE_PROVENANCE_COLUMNS | {"subpanel_image", "layout_plan"}
REQUIRED_SUBPANEL_DIMENSION_COLUMNS = {
    "figure",
    "panel",
    "subpanel_png",
    "width_px",
    "height_px",
    "width_in",
    "height_in",
}
REQUIRED_PANEL_MAP_COLUMNS = {
    "figure",
    "panel",
    "approved_source",
    "approved_generator",
    "intended_content",
    "regeneration_strategy",
}
REQUIRED_LAYOUT_COLUMNS = {
    "figure",
    "panel",
    "subpanel_png",
    "x_in",
    "y_in",
    "width_in",
    "height_in",
    "sx",
    "sy",
    "x_npc",
    "y_npc",
    "width_npc",
    "height_npc",
    "layout_width_in",
    "layout_height_in",
}


class Reporter:
    def __init__(self) -> None:
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.ok: List[str] = []

    def error(self, msg: str) -> None:
        self.errors.append(msg)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)

    def pass_(self, msg: str) -> None:
        self.ok.append(msg)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "status": "ERROR" if self.errors else "OK",
            "errors": self.errors,
            "warnings": self.warnings,
            "ok": self.ok,
        }


def load_contract(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve(root: Path, value: Optional[str]) -> Optional[Path]:
    if not value:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def as_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        return [str(item) for item in value]
    return []


def png_dimensions(path: Path) -> Tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or header[:8] != PNG_SIGNATURE:
        raise ValueError("not a PNG file")
    width, height = struct.unpack(">II", header[16:24])
    return width, height


def check_file_exists(reporter: Reporter, path: Path, label: str) -> bool:
    if not path.exists():
        reporter.error(f"Missing {label}: {path}")
        return False
    if not path.is_file():
        reporter.error(f"{label} is not a file: {path}")
        return False
    if path.stat().st_size <= 0:
        reporter.error(f"{label} is empty: {path}")
        return False
    reporter.pass_(f"Found {label}: {path}")
    return True


def check_png(reporter: Reporter, path: Path, label: str, warn_on_name: bool = False) -> None:
    if not check_file_exists(reporter, path, label):
        return
    try:
        width, height = png_dimensions(path)
    except Exception as exc:  # noqa: BLE001
        reporter.error(f"{label} is not a readable PNG: {path} ({exc})")
        return
    if width <= 0 or height <= 0:
        reporter.error(f"{label} has invalid dimensions {width}x{height}: {path}")
    elif width < 500 or height < 500:
        reporter.warn(f"{label} dimensions are unusually small ({width}x{height}): {path}")
    else:
        reporter.pass_(f"Readable PNG {label}: {path} ({width}x{height})")
    if warn_on_name and re.search(r"(draft|prototype|option)", path.name, re.IGNORECASE):
        reporter.warn(f"Polished output filename still looks provisional: {path.name}")


def check_parent(reporter: Reporter, path: Path, label: str) -> None:
    parent = path if path.suffix == "" else path.parent
    if parent.exists():
        if os.access(parent, os.W_OK):
            reporter.pass_(f"{label} parent is writable: {parent}")
            return
        try:
            with tempfile.NamedTemporaryFile(prefix=".validate_write_", dir=parent):
                pass
        except OSError:
            reporter.error(f"{label} parent is not writable: {parent}")
        else:
            reporter.pass_(f"{label} parent is writable by create/delete check: {parent}")
    elif parent.parent.exists():
        reporter.warn(f"{label} directory does not exist yet but parent exists: {parent}")
    else:
        reporter.error(f"{label} parent directory does not exist: {parent.parent}")


def is_under(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def raster_call_png_literals(script_text: str) -> List[str]:
    literals: List[str] = []
    for match in RASTER_CALL_RE.finditer(script_text):
        literals.extend(PNG_LITERAL_RE.findall(match.group(2)))
    return literals


def discover_approved_raster_root(reporter: Reporter, root: Path) -> Optional[Path]:
    """Find the single centralized approved-raster directory without assuming its parent."""
    candidates: Set[Path] = set()
    if root.name == APPROVED_RASTER_DIRNAME and root.is_dir():
        candidates.add(root.resolve())
    direct = root / APPROVED_RASTER_DIRNAME
    if direct.is_dir():
        candidates.add(direct.resolve())
    try:
        top_level_dirs = [item for item in root.iterdir() if item.is_dir()]
    except OSError as exc:
        reporter.error(f"Could not inspect project root for approved raster directory: {root} ({exc})")
        return None
    for parent in top_level_dirs:
        candidate = parent / APPROVED_RASTER_DIRNAME
        if candidate.is_dir():
            candidates.add(candidate.resolve())

    if not candidates:
        reporter.error(
            "Immutable raster use requires one centralized repository directory named "
            f"{APPROVED_RASTER_DIRNAME!r}; none was found at repository level. "
            "Explain the issue to the user and block raster use."
        )
        return None
    if len(candidates) > 1:
        rendered = "; ".join(str(item) for item in sorted(candidates))
        reporter.error(
            "Immutable raster use requires exactly one centralized repository directory named "
            f"{APPROVED_RASTER_DIRNAME!r}; found {len(candidates)}: {rendered}. "
            "Explain the ambiguity to the user and block raster use."
        )
        return None

    allowed_root = next(iter(candidates))
    reporter.pass_(f"Discovered canonical approved raster directory: {allowed_root}")
    return allowed_root


def check_approved_raster_subpanels(
    reporter: Reporter,
    root: Path,
    contract: Dict[str, Any],
) -> Set[Path]:
    approved: Set[Path] = set()
    declared = as_list(contract.get("approved_raster_subpanels"))
    if not declared:
        return approved
    allowed_root = discover_approved_raster_root(reporter, root)
    for item in declared:
        path = resolve(root, item)
        if path is None:
            continue
        approved.add(path.resolve())
        if allowed_root is None:
            reporter.error(f"Cannot authorize raster subpanel without a canonical directory: {path}")
        elif not is_under(path, allowed_root):
            reporter.error(f"Approved raster subpanel must be under {allowed_root}: {path}")
        check_png(reporter, path, "approved raster subpanel")
    return approved


def check_r_script(reporter: Reporter, path: Path, label: str, phase: str) -> None:
    if path.suffix.lower() != ".r":
        reporter.error(f"{label} must be an R script ending in .R: {path}")
        return
    if phase in {"preflight", "both"}:
        check_parent(reporter, path, label)
    if phase in {"postflight", "both"}:
        check_file_exists(reporter, path, label)


def check_figure_script_content(
    reporter: Reporter,
    root: Path,
    path: Path,
    contract: Dict[str, Any],
    label: str,
    approved_raster_paths: Set[Path],
) -> None:
    try:
        script_text = path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return
    if not re.search(r"\b(cowplot|patchwork|grid|gridExtra)\b", script_text):
        reporter.warn(f"{label} does not mention cowplot, patchwork, grid, or gridExtra")
    if re.search(r"\bgrid\b|viewport\s*\(", script_text) and not re.search(
        r"\b(cowplot|patchwork)\b", script_text
    ):
        reporter.warn(f"{label} appears to use raw grid layout; verify bottom-left layout coordinates in notes")
    if Y_INVERSION_RE.search(script_text):
        reporter.error(
            f"{label} appears to invert layout y coordinates; "
            "layout_plan y_npc is already a lower-left origin"
        )
    if RASTER_ASSEMBLY_RE.search(script_text):
        raster_literals = raster_call_png_literals(script_text)
        if not approved_raster_paths:
            reporter.error(
                f"{label} appears to assemble figures from raster image files, "
                "but the contract declares no approved_raster_subpanels."
            )
        elif not raster_literals:
            reporter.error(
                f"{label} uses raster-loading calls, but their PNG inputs could not be "
                "validated statically. Reference approved_raster_subpanels paths directly "
                "in raster-loading calls."
            )
        else:
            unapproved: List[str] = []
            approved_basenames = {item.name for item in approved_raster_paths}
            for item in raster_literals:
                literal_path = resolve(root, item)
                literal_uses_basename_only = "/" not in item and "\\" not in item
                declared_path = literal_path is not None and literal_path.resolve() in approved_raster_paths
                declared_basename = literal_uses_basename_only and Path(item).name in approved_basenames
                if not declared_path and not declared_basename:
                    unapproved.append(item)
            if unapproved:
                preview = "; ".join(unapproved[:8])
                more = "" if len(unapproved) <= 8 else f"; ... {len(unapproved) - 8} more"
                reporter.error(
                    f"{label} raster-loading calls reference PNG paths not declared "
                    f"in approved_raster_subpanels: {preview}{more}"
                )
            else:
                reporter.pass_(f"{label} raster-loading calls use declared approved raster subpanels")
    if not contract.get("allow_in_figure_titles", False) and IN_FIGURE_TITLE_RE.search(script_text):
        reporter.error(
            f"{label} appears to add figure-level title/subtitle/caption text. "
            "Final PNGs should not contain figure titles, subtitles, captions, or Figure headers."
        )
    if contract.get("layout_plan"):
        layout_basename = Path(str(contract["layout_plan"])).name
        if layout_basename not in script_text:
            reporter.warn(f"{label} does not appear to reference the layout plan")


def read_csv_rows(reporter: Reporter, path: Path, label: str) -> Tuple[List[Dict[str, str]], Set[str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception as exc:  # noqa: BLE001
        reporter.error(f"Could not read {label} CSV: {path} ({exc})")
        return [], set()
    if not rows:
        reporter.error(f"{label} has no rows: {path}")
        return [], set()
    return rows, set(rows[0].keys())


def parse_positive_float(value: str, label: str, reporter: Reporter) -> Optional[float]:
    try:
        parsed = float(value)
    except Exception:  # noqa: BLE001
        reporter.error(f"{label} is not numeric: {value!r}")
        return None
    if parsed <= 0:
        reporter.error(f"{label} must be positive: {value!r}")
        return None
    return parsed


def panel_key(row: Dict[str, Any]) -> Tuple[str, str]:
    return str(row.get("figure", "")).strip(), str(row.get("panel", "")).strip()


def check_paths_under_polish_root(
    reporter: Reporter,
    root: Path,
    polish_root: Optional[Path],
    contract: Dict[str, Any],
) -> None:
    if polish_root is None or not contract.get("require_polish_root", True):
        return
    output_fields = [
        "figure_script",
        "subpanel_script",
        "panel_map",
        "subpanel_dimensions",
        "layout_plan",
        "layout_report",
        "assembly_script",
        "output_manifest",
        "rebuild_manifest",
        "byte_identity_report",
        "polishing_notes",
        "visual_qc_file",
        "validation_report",
        "provenance_table",
        "report_dir",
        "output_dir",
    ]
    output_list_fields = ["expected_outputs", "legend_files", "layout_qc_files", "optional_outputs"]
    for field in output_fields:
        value = contract.get(field)
        if not value:
            continue
        path = resolve(root, str(value))
        if path is not None and not is_under(path, polish_root):
            reporter.error(f"Contract output field {field} is outside polish_root: {path}")
    for field in output_list_fields:
        for item in as_list(contract.get(field)):
            path = resolve(root, item)
            if path is not None and not is_under(path, polish_root):
                reporter.error(f"Contract output field {field} is outside polish_root: {path}")

    final_dir = polish_root / "final_images"
    for item in as_list(contract.get("expected_outputs")):
        path = resolve(root, item)
        if path is not None and not is_under(path, final_dir):
            reporter.warn(f"Polished output is under polish_root but not final_images/: {path}")

    for field, expected_dir in [
        ("figure_script", "scripts"),
        ("subpanel_script", "scripts"),
        ("assembly_script", "scripts"),
        ("subpanel_dimensions", "layout"),
        ("layout_plan", "layout"),
        ("layout_report", "layout"),
    ]:
        value = contract.get(field)
        if not value:
            continue
        path = resolve(root, str(value))
        expected_parent = polish_root / expected_dir
        if path is not None and not is_under(path, expected_parent):
            reporter.warn(f"{field} is not under polish_root/{expected_dir}/: {path}")


def check_feedback_context(reporter: Reporter, root: Path, contract: Dict[str, Any], phase: str) -> None:
    context_value = contract.get("feedback_manager_context")
    if not context_value:
        return
    context_path = resolve(root, str(context_value))
    if context_path is None:
        return
    if phase in {"preflight", "both"}:
        check_parent(reporter, context_path, "feedback-manager context")
    if phase in {"postflight", "both"}:
        check_file_exists(reporter, context_path, "feedback-manager context")


def split_manifest_cell(value: str) -> List[str]:
    parts: List[str] = []
    for chunk in str(value).split(";"):
        chunk = chunk.strip().strip('"').strip("'")
        if chunk:
            parts.append(chunk)
    return parts


def check_manifest_paths(reporter: Reporter, root: Path, manifest: Path) -> None:
    if manifest.suffix.lower() != ".csv":
        reporter.warn(f"Manifest path checking only supports CSV: {manifest}")
        return
    try:
        with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception as exc:  # noqa: BLE001
        reporter.error(f"Could not read CSV manifest: {manifest} ({exc})")
        return
    if not rows:
        reporter.warn(f"CSV manifest has no rows: {manifest}")
        return
    missing: List[str] = []
    checked = 0
    for row in rows:
        for value in row.values():
            for candidate in split_manifest_cell(str(value)):
                if not PATH_LIKE_RE.match(candidate):
                    continue
                path = resolve(root, candidate)
                if path is None:
                    continue
                checked += 1
                if not path.exists():
                    missing.append(candidate)
    if missing:
        preview = "; ".join(missing[:8])
        more = "" if len(missing) <= 8 else f"; ... {len(missing) - 8} more"
        reporter.error(f"Manifest references missing paths: {preview}{more}")
    elif checked:
        reporter.pass_(f"Manifest path references exist ({checked} checked): {manifest}")
    else:
        reporter.warn(f"No local path-like cells found in manifest: {manifest}")


def check_legend_file(reporter: Reporter, path: Path) -> None:
    if not check_file_exists(reporter, path, "legend file"):
        return
    text = path.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        reporter.error(f"Legend file is blank: {path}")
        return
    if not re.search(r"\bFigure\s+S?\d+\s*,", text):
        reporter.warn(f"Legend does not appear to begin with 'Figure N,' or 'Figure SN,': {path}")
    if not re.search(r"\b[a-z]\b", text):
        reporter.warn(f"Legend may be missing lowercase panel-label descriptions: {path}")
    reporter.pass_(f"Legend file is non-empty: {path}")


def check_visual_qc_file(reporter: Reporter, path: Path) -> None:
    if not check_file_exists(reporter, path, "visual QC file"):
        return
    text = path.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        reporter.error(f"Visual QC file is blank: {path}")
        return
    lower = text.lower()
    required_topics = {
        "title/subtitle": ["title", "subtitle", "caption", "header"],
        "clipping": ["clip", "clipping", "cropped", "cut off"],
        "spacing/layout": ["spacing", "layout", "margin", "gutter", "alignment", "overlap"],
        "readability": ["readability", "readable", "text", "font"],
    }
    missing_topics = [
        label for label, terms in required_topics.items() if not any(term in lower for term in terms)
    ]
    if missing_topics:
        reporter.warn(f"Visual QC file may be missing checklist topics: {', '.join(missing_topics)}")
    if not re.search(r"\b(pass|passed|ok|fail|failed|rerun|revised|accepted caveat)\b", lower):
        reporter.warn(f"Visual QC file does not appear to record pass/fail or revision status: {path}")
    reporter.pass_(f"Visual QC file is non-empty: {path}")


def expected_letter_sequence(n: int) -> List[str]:
    return [chr(ord("a") + i) for i in range(n)]


def check_panel_map(reporter: Reporter, root: Path, path: Path) -> Set[Tuple[str, str]]:
    keys: Set[Tuple[str, str]] = set()
    if not check_file_exists(reporter, path, "panel map"):
        return keys
    if path.suffix.lower() != ".csv":
        reporter.error(f"Panel map must be CSV: {path}")
        return keys
    rows, columns = read_csv_rows(reporter, path, "panel map")
    if not rows:
        return keys
    missing_cols = sorted(REQUIRED_PANEL_MAP_COLUMNS - columns)
    if missing_cols:
        reporter.error(f"Panel map missing required columns: {', '.join(missing_cols)}")
        return keys

    duplicate_keys: List[Tuple[str, str]] = []
    empty_required: List[str] = []
    missing_paths: List[str] = []
    panels_by_figure: Dict[str, List[str]] = {}
    for i, row in enumerate(rows, start=2):
        key = panel_key(row)
        row_id = f"row {i} ({key[0] or '?'} {key[1] or '?'})"
        if not key[0] or not key[1]:
            reporter.error(f"Panel map has blank figure/panel in {row_id}")
            continue
        if not re.fullmatch(r"[a-z]", key[1]):
            reporter.error(f"{row_id}: panel label must be one lowercase letter")
        if key in keys:
            duplicate_keys.append(key)
        keys.add(key)
        panels_by_figure.setdefault(key[0], []).append(key[1])

        for col in REQUIRED_PANEL_MAP_COLUMNS:
            if not str(row.get(col, "")).strip():
                empty_required.append(f"{row_id}: {col}")
        for col in ["approved_source", "approved_generator"]:
            for candidate in split_manifest_cell(str(row.get(col, ""))):
                if not PATH_LIKE_RE.match(candidate):
                    continue
                resolved = resolve(root, candidate)
                if resolved is not None and not resolved.exists():
                    missing_paths.append(f"{row_id}: {candidate}")

    if duplicate_keys:
        reporter.error(f"Panel map has duplicate figure/panel keys: {duplicate_keys[:8]}")
    if empty_required:
        preview = "; ".join(empty_required[:8])
        more = "" if len(empty_required) <= 8 else f"; ... {len(empty_required) - 8} more"
        reporter.error(f"Panel map has empty required cells: {preview}{more}")
    if missing_paths:
        preview = "; ".join(missing_paths[:8])
        more = "" if len(missing_paths) <= 8 else f"; ... {len(missing_paths) - 8} more"
        reporter.error(f"Panel map references missing paths: {preview}{more}")

    for figure, panels in panels_by_figure.items():
        expected = expected_letter_sequence(len(panels))
        if panels != expected:
            reporter.warn(
                f"Panel map labels for {figure} are not sequential in row order: {panels}; expected {expected}"
            )

    if keys:
        reporter.pass_(f"Panel map defines {len(keys)} expected displayed panel(s): {path}")
    return keys


def check_panel_key_match(
    reporter: Reporter,
    observed_keys: Set[Tuple[str, str]],
    expected_keys: Set[Tuple[str, str]],
    observed_label: str,
) -> None:
    missing = sorted(expected_keys - observed_keys)
    extra = sorted(observed_keys - expected_keys)
    if missing:
        reporter.error(f"{observed_label} is missing panel-map keys: {missing[:8]}")
    if extra:
        reporter.error(f"{observed_label} has panel keys not present in panel map: {extra[:8]}")
    if expected_keys and not missing and not extra:
        reporter.pass_(f"{observed_label} matches panel map exactly ({len(expected_keys)} panels)")


def check_layout_qc_file(reporter: Reporter, path: Path) -> None:
    if path.suffix.lower() == ".png":
        check_png(reporter, path, "layout QC PNG")
    else:
        check_file_exists(reporter, path, "layout QC file")


def check_subpanel_dimensions(
    reporter: Reporter,
    root: Path,
    path: Path,
    polish_root: Optional[Path] = None,
    expected_keys: Optional[Set[Tuple[str, str]]] = None,
) -> Set[Tuple[str, str]]:
    keys: Set[Tuple[str, str]] = set()
    if not check_file_exists(reporter, path, "subpanel dimension table"):
        return keys
    if path.suffix.lower() != ".csv":
        reporter.error(f"Subpanel dimension table must be CSV: {path}")
        return keys
    rows, columns = read_csv_rows(reporter, path, "subpanel dimension table")
    if not rows:
        return keys
    missing_cols = sorted(REQUIRED_SUBPANEL_DIMENSION_COLUMNS - columns)
    if missing_cols:
        reporter.error(f"Subpanel dimension table missing required columns: {', '.join(missing_cols)}")
        return keys

    duplicate_keys: List[Tuple[str, str]] = []
    for i, row in enumerate(rows, start=2):
        key = panel_key(row)
        row_id = f"row {i} ({key[0] or '?'} {key[1] or '?'})"
        if not key[0] or not key[1]:
            reporter.error(f"Subpanel dimension table has blank figure/panel in {row_id}")
            continue
        if key in keys:
            duplicate_keys.append(key)
        keys.add(key)
        width_px = parse_positive_float(str(row.get("width_px", "")), f"{row_id}: width_px", reporter)
        height_px = parse_positive_float(str(row.get("height_px", "")), f"{row_id}: height_px", reporter)
        parse_positive_float(str(row.get("width_in", "")), f"{row_id}: width_in", reporter)
        parse_positive_float(str(row.get("height_in", "")), f"{row_id}: height_in", reporter)
        png_path = resolve(root, str(row.get("subpanel_png", "")).strip())
        if png_path is None:
            reporter.error(f"{row_id}: subpanel_png is blank")
            continue
        check_png(reporter, png_path, f"regenerated subpanel PNG {key[0]} {key[1]}")
        if polish_root is not None and not is_under(png_path, polish_root / "subpanels"):
            reporter.warn(f"{row_id}: regenerated subpanel is not under polish_root/subpanels/: {png_path}")
        if png_path.exists():
            try:
                actual_width, actual_height = png_dimensions(png_path)
                if width_px is not None and int(round(width_px)) != actual_width:
                    reporter.error(
                        f"{row_id}: width_px {width_px:g} does not match PNG width {actual_width}: {png_path}"
                    )
                if height_px is not None and int(round(height_px)) != actual_height:
                    reporter.error(
                        f"{row_id}: height_px {height_px:g} does not match PNG height {actual_height}: {png_path}"
                    )
            except Exception:
                pass
    if duplicate_keys:
        reporter.error(f"Subpanel dimension table has duplicate figure/panel keys: {duplicate_keys[:8]}")
    if expected_keys is not None:
        check_panel_key_match(reporter, keys, expected_keys, "Subpanel dimension table")
    if keys:
        reporter.pass_(f"Subpanel dimension table covers {len(keys)} regenerated subpanel(s): {path}")
    return keys


def check_layout_plan(
    reporter: Reporter,
    root: Path,
    path: Path,
    expected_keys: Optional[Set[Tuple[str, str]]] = None,
) -> Set[Tuple[str, str]]:
    keys: Set[Tuple[str, str]] = set()
    if not check_file_exists(reporter, path, "layout plan"):
        return keys
    if path.suffix.lower() != ".csv":
        reporter.error(f"Layout plan must be CSV: {path}")
        return keys
    rows, columns = read_csv_rows(reporter, path, "layout plan")
    if not rows:
        return keys
    missing_cols = sorted(REQUIRED_LAYOUT_COLUMNS - columns)
    if missing_cols:
        reporter.error(f"Layout plan missing required columns: {', '.join(missing_cols)}")
        return keys

    positive_cols = [
        "width_in",
        "height_in",
        "sx",
        "sy",
        "width_npc",
        "height_npc",
        "layout_width_in",
        "layout_height_in",
    ]
    nonnegative_cols = ["x_in", "y_in", "x_npc", "y_npc"]
    duplicate_keys: List[Tuple[str, str]] = []
    for i, row in enumerate(rows, start=2):
        key = panel_key(row)
        row_id = f"row {i} ({key[0] or '?'} {key[1] or '?'})"
        if not key[0] or not key[1]:
            reporter.error(f"Layout plan has blank figure/panel in {row_id}")
            continue
        if key in keys:
            duplicate_keys.append(key)
        keys.add(key)
        for col in positive_cols:
            parse_positive_float(str(row.get(col, "")), f"{row_id}: {col}", reporter)
        for col in nonnegative_cols:
            try:
                value = float(str(row.get(col, "")))
            except Exception:  # noqa: BLE001
                reporter.error(f"{row_id}: {col} is not numeric: {row.get(col, '')!r}")
                continue
            if value < 0:
                reporter.error(f"{row_id}: {col} must be nonnegative: {value:g}")
        subpanel = resolve(root, str(row.get("subpanel_png", "")).strip())
        if subpanel is None:
            reporter.error(f"{row_id}: subpanel_png is blank")
        elif not subpanel.exists():
            reporter.error(f"{row_id}: layout plan references missing subpanel_png: {subpanel}")

    if duplicate_keys:
        reporter.error(f"Layout plan has duplicate figure/panel keys: {duplicate_keys[:8]}")
    if expected_keys is not None:
        missing = sorted(expected_keys - keys)
        extra = sorted(keys - expected_keys)
        if missing:
            reporter.error(f"Layout plan is missing regenerated subpanel keys: {missing[:8]}")
        if extra:
            reporter.error(f"Layout plan has panel keys not present in expected panel set: {extra[:8]}")
        if not missing and not extra:
            reporter.pass_(f"Layout plan matches expected panel set exactly ({len(expected_keys)} panels)")
    if keys:
        reporter.pass_(f"Layout plan covers {len(keys)} panel placement(s): {path}")
    return keys


def check_provenance_table(
    reporter: Reporter,
    root: Path,
    path: Path,
    require_regenerated: bool,
    expected_keys: Optional[Set[Tuple[str, str]]] = None,
) -> None:
    if not check_file_exists(reporter, path, "provenance table"):
        return
    if path.suffix.lower() != ".csv":
        reporter.error(f"Provenance table must be CSV: {path}")
        return
    rows, columns = read_csv_rows(reporter, path, "provenance table")
    if not rows:
        return
    required_columns = REGENERATED_PROVENANCE_COLUMNS if require_regenerated else LEGACY_PROVENANCE_COLUMNS
    missing_cols = sorted(required_columns - columns)
    if missing_cols:
        reporter.error(f"Provenance table missing required columns: {', '.join(missing_cols)}")
        return

    empty_required: List[str] = []
    missing_paths: List[str] = []
    keys: Set[Tuple[str, str]] = set()
    duplicate_keys: List[Tuple[str, str]] = []
    for i, row in enumerate(rows, start=2):
        row_id = f"row {i} ({row.get('figure', '?')} {row.get('panel', '?')})"
        key = panel_key(row)
        if key in keys:
            duplicate_keys.append(key)
        keys.add(key)
        for col in required_columns:
            if not str(row.get(col, "")).strip():
                empty_required.append(f"{row_id}: {col}")
        for col in ["source_image", "subpanel_image", "output_image", "data_inputs", "layout_plan"]:
            if col not in row:
                continue
            for candidate in split_manifest_cell(str(row.get(col, ""))):
                if not PATH_LIKE_RE.match(candidate):
                    continue
                resolved = resolve(root, candidate)
                if resolved is not None and not resolved.exists():
                    missing_paths.append(f"{row_id}: {candidate}")

    if empty_required:
        preview = "; ".join(empty_required[:8])
        more = "" if len(empty_required) <= 8 else f"; ... {len(empty_required) - 8} more"
        reporter.error(f"Provenance table has empty required cells: {preview}{more}")
    if missing_paths:
        preview = "; ".join(missing_paths[:8])
        more = "" if len(missing_paths) <= 8 else f"; ... {len(missing_paths) - 8} more"
        reporter.error(f"Provenance table references missing paths: {preview}{more}")
    if duplicate_keys:
        reporter.error(f"Provenance table has duplicate figure/panel keys: {duplicate_keys[:8]}")
    if expected_keys is not None:
        missing_keys = sorted(expected_keys - keys)
        extra_keys = sorted(keys - expected_keys)
        if missing_keys:
            reporter.error(f"Provenance table is missing regenerated subpanel keys: {missing_keys[:8]}")
        if extra_keys:
            reporter.error(f"Provenance table has panel keys not present in expected panel set: {extra_keys[:8]}")
        if not missing_keys and not extra_keys:
            reporter.pass_(f"Provenance table matches expected panel set exactly ({len(expected_keys)} panels)")
    if not empty_required and not missing_paths:
        reporter.pass_(f"Provenance table is complete ({len(rows)} rows): {path}")


def run_validation(contract: Dict[str, Any], phase: str) -> Reporter:
    reporter = Reporter()
    root = resolve(Path.cwd(), contract.get("project_root")) or Path.cwd()
    root = root.resolve()
    require_regenerated = contract.get("require_regenerated_subpanels", True)
    require_panel_contract = contract.get("require_panel_contract", require_regenerated)
    require_layout_qc = contract.get("require_layout_qc", True)
    require_visual_qc = contract.get("require_visual_qc", True)
    allow_split_scripts = contract.get("allow_split_scripts", False)

    if "allow_raster_assembly" in contract:
        reporter.error(
            "Contract field allow_raster_assembly has been removed; use "
            "approved_raster_subpanels for raster inputs allowed by the global policy."
        )

    if not contract.get("wp_id"):
        reporter.error("Contract missing required field: wp_id")
    else:
        reporter.pass_(f"WP id: {contract['wp_id']}")

    polish_root = resolve(root, contract.get("polish_root"))
    if contract.get("require_polish_root", True) and polish_root is None:
        reporter.error("Contract missing required field: polish_root")
    elif polish_root is not None:
        if phase in {"preflight", "both"}:
            check_parent(reporter, polish_root, "polish root")
        if phase in {"postflight", "both"}:
            if not polish_root.exists():
                reporter.error(f"Missing polish root: {polish_root}")
            elif not polish_root.is_dir():
                reporter.error(f"polish_root is not a directory: {polish_root}")
            else:
                reporter.pass_(f"Found polish root: {polish_root}")
    check_paths_under_polish_root(reporter, root, polish_root, contract)

    if contract.get("feedback_file"):
        reporter.error(
            "Contract uses deprecated feedback_file; record local revision context in notes "
            "and use feedback_manager_context only for a supplied managed-feedback handoff"
        )
    check_feedback_context(reporter, root, contract, phase)

    panel_map_keys: Set[Tuple[str, str]] = set()
    panel_map = resolve(root, contract.get("panel_map"))
    if require_panel_contract and panel_map is None:
        reporter.error("Contract missing required field: panel_map")
    elif panel_map is not None:
        panel_map_keys = check_panel_map(reporter, root, panel_map)

    for field, label in [
        ("draft_doc", "draft documentation"),
        ("source_files", "source file"),
    ]:
        for item in as_list(contract.get(field)):
            path = resolve(root, item)
            if path is not None:
                check_file_exists(reporter, path, label)

    for item in as_list(contract.get("source_pngs")):
        path = resolve(root, item)
        if path is not None:
            check_png(reporter, path, "source PNG")

    approved_raster_paths = check_approved_raster_subpanels(reporter, root, contract)

    default_optimizer = str(Path(__file__).resolve().parent / "optimize_panel_layout.R")
    layout_optimizer_script = resolve(root, contract.get("layout_optimizer_script") or default_optimizer)
    if require_regenerated and layout_optimizer_script is not None:
        if layout_optimizer_script.suffix.lower() != ".r":
            reporter.error(f"layout optimizer script must be an R script ending in .R: {layout_optimizer_script}")
        else:
            check_file_exists(reporter, layout_optimizer_script, "layout optimizer script")

    if phase in {"preflight", "both"}:
        if require_panel_contract and not contract.get("panel_map"):
            reporter.error("Contract missing required preflight field: panel_map")
        if not allow_split_scripts and (contract.get("subpanel_script") or contract.get("assembly_script")):
            reporter.error(
                "Contract uses split subpanel_script/assembly_script fields; "
                "use one figure_script unless allow_split_scripts is true"
            )
        if require_regenerated:
            if allow_split_scripts:
                for field in ["subpanel_script", "assembly_script"]:
                    if not contract.get(field):
                        reporter.error(f"Contract missing required preflight field: {field}")
            elif not contract.get("figure_script"):
                reporter.error("Contract missing required preflight field: figure_script")
            for field in ["subpanel_dimensions", "layout_plan"]:
                if not contract.get(field):
                    reporter.error(f"Contract missing required preflight field: {field}")
            if not str(contract.get("layout_optimizer_command", "")).strip():
                reporter.error("Contract missing required preflight field: layout_optimizer_command")
        if require_layout_qc and not as_list(contract.get("layout_qc_files")):
            reporter.error("Contract missing required preflight field: layout_qc_files")
        if require_visual_qc and not contract.get("visual_qc_file"):
            reporter.error("Contract missing required preflight field: visual_qc_file")
        for field, label in [
            ("output_dir", "output directory"),
            ("report_dir", "report directory"),
            ("output_manifest", "manifest"),
            ("polishing_notes", "polishing notes"),
            ("visual_qc_file", "visual QC file"),
            ("validation_report", "validation report"),
            ("provenance_table", "provenance table"),
            ("panel_map", "panel map"),
            ("subpanel_dimensions", "subpanel dimension table"),
            ("layout_plan", "layout plan"),
            ("layout_report", "layout report"),
        ]:
            path = resolve(root, contract.get(field))
            if path is not None:
                check_parent(reporter, path, label)
        for field, label in [
            ("figure_script", "figure polishing script"),
            ("subpanel_script", "subpanel generation script"),
            ("assembly_script", "figure assembly script"),
        ]:
            path = resolve(root, contract.get(field))
            required_script = (
                require_regenerated
                and (
                    (field == "figure_script" and not allow_split_scripts)
                    or (field in {"subpanel_script", "assembly_script"} and allow_split_scripts)
                )
            )
            if path is None and required_script:
                reporter.error(f"Contract missing required preflight field: {field}")
            elif path is not None:
                check_r_script(reporter, path, label, "preflight")
        for item in as_list(contract.get("legend_files")):
            path = resolve(root, item)
            if path is not None:
                check_parent(reporter, path, "legend file")
        for item in as_list(contract.get("layout_qc_files")):
            path = resolve(root, item)
            if path is not None:
                check_parent(reporter, path, "layout QC file")
        for item in as_list(contract.get("expected_outputs")):
            path = resolve(root, item)
            if path is not None:
                check_parent(reporter, path, "expected output")

    if phase in {"postflight", "both"}:
        expected = as_list(contract.get("expected_outputs"))
        if not expected:
            reporter.error("Contract missing required postflight field: expected_outputs")
        for item in expected:
            path = resolve(root, item)
            if path is not None:
                check_png(reporter, path, "polished output", warn_on_name=True)

        subpanel_keys: Set[Tuple[str, str]] = set()
        if require_regenerated:
            if not allow_split_scripts and (contract.get("subpanel_script") or contract.get("assembly_script")):
                reporter.error(
                    "Contract uses split subpanel_script/assembly_script fields; "
                    "use one figure_script unless allow_split_scripts is true"
                )
            required_fields = (
                ["subpanel_script", "subpanel_dimensions", "layout_plan", "assembly_script"]
                if allow_split_scripts
                else ["figure_script", "subpanel_dimensions", "layout_plan"]
            )
            for field in required_fields:
                if not contract.get(field):
                    reporter.error(f"Contract missing required postflight field: {field}")

            figure_script = resolve(root, contract.get("figure_script"))
            if figure_script is not None:
                check_r_script(reporter, figure_script, "figure polishing script", "postflight")
                check_figure_script_content(
                    reporter,
                    root,
                    figure_script,
                    contract,
                    "Figure polishing script",
                    approved_raster_paths,
                )

            if allow_split_scripts:
                subpanel_script = resolve(root, contract.get("subpanel_script"))
                if subpanel_script is not None:
                    check_r_script(reporter, subpanel_script, "subpanel generation script", "postflight")

                assembly_script = resolve(root, contract.get("assembly_script"))
                if assembly_script is not None:
                    check_r_script(reporter, assembly_script, "figure assembly script", "postflight")
                    check_figure_script_content(
                        reporter,
                        root,
                        assembly_script,
                        contract,
                        "Figure assembly script",
                        approved_raster_paths,
                    )

            if not str(contract.get("layout_optimizer_command", "")).strip():
                reporter.error("Contract missing required postflight field: layout_optimizer_command")
            else:
                reporter.pass_("Layout optimizer command recorded in contract")

            dimensions = resolve(root, contract.get("subpanel_dimensions"))
            if dimensions is not None:
                subpanel_keys = check_subpanel_dimensions(
                    reporter,
                    root,
                    dimensions,
                    polish_root,
                    panel_map_keys if require_panel_contract and panel_map_keys else None,
                )

            layout_plan = resolve(root, contract.get("layout_plan"))
            if layout_plan is not None:
                expected_layout_keys = panel_map_keys if require_panel_contract and panel_map_keys else subpanel_keys
                check_layout_plan(reporter, root, layout_plan, expected_layout_keys)

            layout_report = resolve(root, contract.get("layout_report"))
            if layout_report is not None:
                check_file_exists(reporter, layout_report, "layout report")
        else:
            subpanel_keys = set()

        legend_files = as_list(contract.get("legend_files"))
        if contract.get("require_legend_files", False) and not legend_files:
            reporter.error("Contract missing required postflight field: legend_files")
        for item in legend_files:
            path = resolve(root, item)
            if path is not None:
                check_legend_file(reporter, path)

        layout_qc_files = as_list(contract.get("layout_qc_files"))
        if require_layout_qc and not layout_qc_files:
            reporter.error("Contract missing required postflight field: layout_qc_files")
        for item in layout_qc_files:
            path = resolve(root, item)
            if path is not None:
                check_layout_qc_file(reporter, path)

        provenance = resolve(root, contract.get("provenance_table"))
        if contract.get("require_provenance_table", True) and provenance is None:
            reporter.error("Contract missing required postflight field: provenance_table")
        elif provenance is not None:
            expected_keys = (
                panel_map_keys
                if require_panel_contract and panel_map_keys
                else subpanel_keys
                if require_regenerated and subpanel_keys
                else None
            )
            check_provenance_table(reporter, root, provenance, require_regenerated, expected_keys)

        manifest = resolve(root, contract.get("output_manifest"))
        if manifest is not None:
            if check_file_exists(reporter, manifest, "output manifest"):
                if contract.get("check_manifest_paths", True):
                    check_manifest_paths(reporter, root, manifest)
        else:
            reporter.warn("No output_manifest specified")

        notes = resolve(root, contract.get("polishing_notes"))
        if notes is not None:
            if check_file_exists(reporter, notes, "polishing notes"):
                text = notes.read_text(encoding="utf-8", errors="replace").lower()
                if contract.get("require_project_map_decision", True):
                    if "project_map" not in text and "project-map" not in text:
                        reporter.warn("Polishing notes do not mention a project-map decision")
                    else:
                        reporter.pass_("Polishing notes mention project-map decision")
        else:
            reporter.warn("No polishing_notes specified")

        visual_qc = resolve(root, contract.get("visual_qc_file"))
        if require_visual_qc and visual_qc is None:
            reporter.error("Contract missing required postflight field: visual_qc_file")
        elif visual_qc is not None:
            check_visual_qc_file(reporter, visual_qc)

        for item in as_list(contract.get("optional_outputs")):
            path = resolve(root, item)
            if path is not None:
                if path.exists():
                    reporter.pass_(f"Optional output present: {path}")
                else:
                    reporter.warn(f"Optional output missing: {path}")

    return reporter


def print_human(report: Dict[str, Any]) -> None:
    for key, label in [("errors", "ERROR"), ("warnings", "WARN"), ("ok", "OK")]:
        for item in report[key]:
            print(f"{label}: {item}")
    print(f"STATUS: {report['status']}")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("contract", type=Path, help="Figure polishing contract JSON")
    parser.add_argument(
        "--phase",
        choices=["preflight", "postflight", "both"],
        default="both",
        help="Validation phase to run",
    )
    parser.add_argument("--json", action="store_true", help="Print JSON report")
    parser.add_argument("--write-report", type=Path, help="Write JSON report to this path")
    args = parser.parse_args(argv)

    try:
        contract = load_contract(args.contract)
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: Could not load contract: {args.contract} ({exc})", file=sys.stderr)
        return 1

    reporter = run_validation(contract, args.phase)
    report = reporter.as_dict()
    report["contract"] = str(args.contract)
    report["phase"] = args.phase

    if args.write_report:
        args.write_report.parent.mkdir(parents=True, exist_ok=True)
        args.write_report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print_human(report)

    return 1 if reporter.errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
