#!/usr/bin/env python
"""Batch-convert an entire pypgx-bundle ``cnv/`` tree of pickled ``Model[CNV]``
archives to the Rust-native (``data.json``) form, and optionally verify that each
converted model reproduces scikit-learn's ``predict`` exactly.

PyPGx ships its CNV callers as pickled sklearn ``OneVsRestClassifier(SVC)``
objects, which cannot be loaded in Rust. This driver runs the same extraction as
``convert_cnv_model.py`` (support vectors, dual coefs, intercept, ``_gamma``) over
every ``{assembly}/{gene}.zip`` under a source ``cnv/`` directory, writing a
mirror tree of Rust-readable archives. The Rust ``pypgx::api::predict_cnv`` then
reproduces sklearn's predictions natively (see ``src/cnv.rs`` / ``tests/cnv.rs``).

Parity target: sklearn *in this environment* is exactly what PyPGx itself uses to
load these same pickles, so ``--verify`` compares the extracted-param decision
against ``clf.predict`` over random copy-number vectors — 0 mismatches means the
weights are faithful.

Usage:
    python tools/convert_cnv_models_all.py SRC_cnv_dir DST_cnv_dir [--verify [N]]

Run inside ``.refenv`` (or any env with pypgx + scikit-learn). Example:
    python tools/convert_cnv_models_all.py ~/pypgx-bundle/cnv ./rust-bundle/cnv --verify
"""
import sys
import os
import json
import glob
import zipfile
import pickle
import warnings

warnings.filterwarnings("ignore")  # cross-version unpickle warning is expected


def load_pickled_model(path: str):
    """Return (clf, metadata_text) from a PyPGx pickled Model[CNV] archive.

    No pypgx import needed: the archive is a zip of `<stem>/metadata.txt` +
    `<stem>/data.sav`, and `data.sav` is a pure sklearn `OneVsRestClassifier`
    pickle (only scikit-learn + numpy are needed to unpickle it)."""
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        sav = next(n for n in names if n.endswith("data.sav"))
        meta = next((n for n in names if n.endswith("metadata.txt")), None)
        clf = pickle.loads(z.read(sav))
        metadata = z.read(meta).decode() if meta else ""
    return clf, metadata


def extract(clf) -> dict:
    """Pull each binary RBF-SVM estimator's fitted params into a plain dict that
    matches ``pypgx::cnv::CnvModel`` (serde)."""
    estimators = []
    for est, label in zip(clf.estimators_, clf.classes_):
        estimators.append({
            "label": int(label),
            "gamma": float(est._gamma),
            "intercept": float(est.intercept_[0]),
            "dual_coef": est.dual_coef_[0].tolist(),
            "support_vectors": est.support_vectors_.tolist(),
        })
    return {"classes": clf.classes_.tolist(), "estimators": estimators}


def write_archive(out_path: str, metadata_text: str, model: dict) -> None:
    stem = os.path.splitext(os.path.basename(out_path))[0]
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{stem}/metadata.txt", metadata_text)
        z.writestr(f"{stem}/data.json", json.dumps(model))


def verify(clf, model: dict, n: int) -> int:
    """Compare clf.predict vs a numpy reimplementation from the extracted params
    over ``n`` random copy-number vectors; return the mismatch count."""
    import numpy as np

    rng = np.random.default_rng(0)
    ndim = clf.estimators_[0].support_vectors_.shape[1]
    X = rng.uniform(0.0, 4.0, size=(n, ndim))
    scores = np.empty((n, len(model["estimators"])))
    for k, e in enumerate(model["estimators"]):
        sv = np.asarray(e["support_vectors"])
        dc = np.asarray(e["dual_coef"])
        sq = ((sv[None, :, :] - X[:, None, :]) ** 2).sum(-1)
        scores[:, k] = (dc[None, :] * np.exp(-e["gamma"] * sq)).sum(1) + e["intercept"]
    classes = np.asarray(model["classes"])
    return int((clf.predict(X) != classes[scores.argmax(1)]).sum())


def main() -> None:
    args = sys.argv[1:]
    do_verify = False
    n = 300
    if "--verify" in args:
        i = args.index("--verify")
        args.pop(i)
        if i < len(args) and args[i].isdigit():
            n = int(args.pop(i))
        do_verify = True
    if len(args) != 2:
        sys.exit(__doc__)
    src, dst = args

    total_mismatch = 0
    count = 0
    for path in sorted(glob.glob(f"{src}/*/*.zip")):
        assembly = os.path.basename(os.path.dirname(path))
        gene = os.path.splitext(os.path.basename(path))[0]
        clf, metadata_text = load_pickled_model(path)
        model = extract(clf)
        os.makedirs(f"{dst}/{assembly}", exist_ok=True)
        write_archive(f"{dst}/{assembly}/{gene}.zip", metadata_text, model)
        count += 1
        line = f"{assembly:7} {gene:9} {len(model['classes']):3d} classes"
        if do_verify:
            m = verify(clf, model, n)
            total_mismatch += m
            line += f"  mismatch/{n}: {m}" + ("  <-- MISMATCH" if m else "")
        print(line)

    print(f"\n{count} models converted -> {dst}")
    if do_verify:
        print(f"TOTAL prediction mismatches across all models: {total_mismatch}")
        if total_mismatch:
            sys.exit(1)


if __name__ == "__main__":
    main()
