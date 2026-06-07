#!/usr/bin/env python
"""Convert a PyPGx pickled ``Model[CNV]`` archive to the Rust-native form.

PyPGx ships its CNV callers as pickled scikit-learn ``OneVsRestClassifier(SVC)``
objects (``data.sav`` inside the archive), which cannot be loaded in Rust. This
script extracts each binary estimator's fitted RBF-SVM parameters (support
vectors, dual coefficients, intercept, gamma) and writes a new archive whose
``data.json`` matches ``pypgx::cnv::CnvModel``. The Rust ``predict_cnv`` then
reproduces sklearn's predictions exactly (verified in ``tests/cnv.rs``).

Usage (in any env with scikit-learn that can unpickle the model):
    python tools/convert_cnv_model.py IN_pickled_Model.zip OUT_rust_Model.zip

No pypgx import is needed — the archive is a zip of `<stem>/metadata.txt` +
`<stem>/data.sav`, and `data.sav` is a pure sklearn `OneVsRestClassifier(SVC)`
pickle, so only scikit-learn (+ numpy) are required to unpickle it. Point it at a
model from the ``pypgx-bundle`` (e.g. ``cnv/GRCh37/CYP2D6.zip``).
"""
import sys
import json
import zipfile
import pickle
import os
import warnings

warnings.filterwarnings("ignore")  # cross-version unpickle warning is expected


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    in_path, out_path = sys.argv[1], sys.argv[2]

    # Read the pickled sklearn model + metadata straight out of the zip.
    with zipfile.ZipFile(in_path) as z:
        names = z.namelist()
        sav = next(n for n in names if n.endswith("data.sav"))
        meta = next((n for n in names if n.endswith("metadata.txt")), None)
        clf = pickle.loads(z.read(sav))  # OneVsRestClassifier(SVC(kernel='rbf'))
        metadata = z.read(meta).decode() if meta else ""

    estimators = []
    for est, label in zip(clf.estimators_, clf.classes_):
        estimators.append({
            "label": int(label),
            "gamma": float(est._gamma),
            "intercept": float(est.intercept_[0]),
            "dual_coef": est.dual_coef_[0].tolist(),
            "support_vectors": est.support_vectors_.tolist(),
        })
    model = {"classes": clf.classes_.tolist(), "estimators": estimators}

    # Write a Rust-readable archive: <stem>/metadata.txt + <stem>/data.json.
    stem = os.path.splitext(os.path.basename(out_path))[0]
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{stem}/metadata.txt", metadata)
        z.writestr(f"{stem}/data.json", json.dumps(model))
    print(f"Wrote Rust Model[CNV] -> {out_path} "
          f"({len(estimators)} estimators, {len(model['classes'])} classes)")


if __name__ == "__main__":
    main()
