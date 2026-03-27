#!/usr/bin/env python3
"""
Train a simple object-level classifier from the embedding feature table.

This first-pass trainer uses:
- target: `label_int`
- features: `log_area_px` + `obj_emb_*` + `frame_emb_*`

It performs a group holdout by frame by default so objects from the same frame
do not appear in both train and test splits.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
from typing import List

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    log_loss,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit, StratifiedShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train a simple object embedding classifier.")
    parser.add_argument("--input_csv", required=True, help="Object feature table CSV.")
    parser.add_argument("--output_dir", required=True, help="Directory for model and evaluation outputs.")
    parser.add_argument("--target_col", default="label_int", help="Binary target column.")
    parser.add_argument("--group_col", default="frame", help="Grouping column for held-out split.")
    parser.add_argument("--test_size", type=float, default=0.2, help="Held-out fraction.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed.")
    parser.add_argument("--max_iter", type=int, default=2000, help="Logistic regression max iterations.")
    parser.add_argument(
        "--no_group_split",
        action="store_true",
        help="Use a row-level stratified split instead of holding out complete frames.",
    )
    return parser.parse_args()


def get_feature_columns(df: pd.DataFrame) -> List[str]:
    feature_cols = ["log_area_px"]
    feature_cols.extend(sorted(c for c in df.columns if c.startswith("obj_emb_")))
    feature_cols.extend(sorted(c for c in df.columns if c.startswith("frame_emb_")))
    missing = [c for c in ["log_area_px"] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required feature columns: {missing}")
    if len(feature_cols) <= 1:
        raise ValueError("No embedding columns found. Expected obj_emb_* and frame_emb_* columns.")
    return feature_cols


def split_indices(
    df: pd.DataFrame,
    target: np.ndarray,
    group_col: str,
    test_size: float,
    seed: int,
    no_group_split: bool,
) -> tuple[np.ndarray, np.ndarray, str]:
    idx = np.arange(df.shape[0], dtype=np.int64)
    if no_group_split:
        splitter = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
        train_idx, test_idx = next(splitter.split(idx.reshape(-1, 1), target))
        return train_idx, test_idx, "row_stratified"

    if group_col not in df.columns:
        raise ValueError(f"group_col '{group_col}' not found in input data")
    groups = df[group_col].to_numpy()
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
    train_idx, test_idx = next(splitter.split(idx, target, groups=groups))
    return train_idx, test_idx, f"group_by_{group_col}"


def metric_dict(y_true: np.ndarray, y_prob: np.ndarray, threshold: float = 0.5) -> dict:
    y_pred = (y_prob >= threshold).astype(np.int32)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    metrics = {
        "n": int(y_true.size),
        "prevalence": float(np.mean(y_true)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "threshold": float(threshold),
    }
    if np.unique(y_true).size > 1:
        metrics["roc_auc"] = float(roc_auc_score(y_true, y_prob))
        metrics["average_precision"] = float(average_precision_score(y_true, y_prob))
        metrics["log_loss"] = float(log_loss(y_true, np.clip(y_prob, 1e-6, 1 - 1e-6)))
    else:
        metrics["roc_auc"] = None
        metrics["average_precision"] = None
        metrics["log_loss"] = None
    return metrics


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    if args.target_col not in df.columns:
        raise ValueError(f"target_col '{args.target_col}' not found in input data")

    feature_cols = get_feature_columns(df)
    use_cols = feature_cols + [args.target_col]
    if not args.no_group_split:
        use_cols.append(args.group_col)

    model_df = df[use_cols].copy()
    model_df = model_df.dropna(axis=0).reset_index(drop=True)

    y = model_df[args.target_col].astype(int).to_numpy()
    X = model_df[feature_cols].astype(np.float32).to_numpy()

    train_idx, test_idx, split_type = split_indices(
        df=model_df,
        target=y,
        group_col=args.group_col,
        test_size=args.test_size,
        seed=args.seed,
        no_group_split=args.no_group_split,
    )

    X_train = X[train_idx]
    X_test = X[test_idx]
    y_train = y[train_idx]
    y_test = y[test_idx]

    model = Pipeline(
        steps=[
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    class_weight="balanced",
                    max_iter=args.max_iter,
                    random_state=args.seed,
                    solver="lbfgs",
                ),
            ),
        ]
    )
    model.fit(X_train, y_train)

    train_prob = model.predict_proba(X_train)[:, 1]
    test_prob = model.predict_proba(X_test)[:, 1]

    train_metrics = metric_dict(y_train, train_prob)
    test_metrics = metric_dict(y_test, test_prob)

    coef = model.named_steps["clf"].coef_.ravel()
    coef_df = pd.DataFrame(
        {
            "feature": feature_cols,
            "coefficient": coef.astype(np.float64, copy=False),
            "abs_coefficient": np.abs(coef).astype(np.float64, copy=False),
        }
    ).sort_values("abs_coefficient", ascending=False, ignore_index=True)
    coef_df.to_csv(os.path.join(args.output_dir, "coefficients.csv"), index=False)

    pred_df = pd.DataFrame(
        {
            "split": ["train"] * len(train_idx) + ["test"] * len(test_idx),
            "row_index": np.concatenate([train_idx, test_idx]).astype(np.int64),
            "y_true": np.concatenate([y_train, y_test]).astype(np.int32),
            "y_prob": np.concatenate([train_prob, test_prob]).astype(np.float64),
        }
    )
    if args.group_col in model_df.columns:
        pred_df["group"] = model_df.iloc[pred_df["row_index"]][args.group_col].to_numpy()
    pred_df["y_pred"] = (pred_df["y_prob"] >= 0.5).astype(np.int32)
    pred_df.to_csv(os.path.join(args.output_dir, "heldout_predictions.csv"), index=False)

    model_bundle = {
        "model": model,
        "feature_columns": feature_cols,
        "target_col": args.target_col,
        "group_col": None if args.no_group_split else args.group_col,
        "split_type": split_type,
        "train_size": int(len(train_idx)),
        "test_size": int(len(test_idx)),
    }
    with open(os.path.join(args.output_dir, "object_embedding_classifier.pkl"), "wb") as handle:
        pickle.dump(model_bundle, handle)

    metrics = {
        "input_csv": os.path.abspath(args.input_csv),
        "output_dir": os.path.abspath(args.output_dir),
        "target_col": args.target_col,
        "group_col": None if args.no_group_split else args.group_col,
        "split_type": split_type,
        "seed": int(args.seed),
        "test_size": float(args.test_size),
        "max_iter": int(args.max_iter),
        "n_rows_input": int(df.shape[0]),
        "n_rows_used": int(model_df.shape[0]),
        "n_features": int(len(feature_cols)),
        "feature_blocks": {
            "log_area_px": 1,
            "obj_emb": int(sum(c.startswith("obj_emb_") for c in feature_cols)),
            "frame_emb": int(sum(c.startswith("frame_emb_") for c in feature_cols)),
        },
        "train_metrics": train_metrics,
        "test_metrics": test_metrics,
    }
    with open(os.path.join(args.output_dir, "metrics.json"), "w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote model: {os.path.abspath(os.path.join(args.output_dir, 'object_embedding_classifier.pkl'))}")
    print(f"Wrote metrics: {os.path.abspath(os.path.join(args.output_dir, 'metrics.json'))}")
    print(f"Wrote predictions: {os.path.abspath(os.path.join(args.output_dir, 'heldout_predictions.csv'))}")
    print(f"Wrote coefficients: {os.path.abspath(os.path.join(args.output_dir, 'coefficients.csv'))}")
    print(
        "Held-out metrics:"
        f" accuracy={test_metrics['accuracy']:.4f}"
        f" balanced_accuracy={test_metrics['balanced_accuracy']:.4f}"
        f" roc_auc={test_metrics['roc_auc'] if test_metrics['roc_auc'] is not None else 'NA'}"
        f" average_precision={test_metrics['average_precision'] if test_metrics['average_precision'] is not None else 'NA'}"
    )


if __name__ == "__main__":
    main()
