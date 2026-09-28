"""
webapp.py — Local Flask web frontend: type a PDB ID in the browser, get
the same annotation report as the CLI, plus an interactive 3D structure
viewer (3Dmol.js), with nothing to install beyond `pip install flask`.

Run with:
    python -m protein_toolkit.webapp
then open http://127.0.0.1:5000 in a browser.

This is meant for local/personal use (no auth, no production WSGI server,
debug reloader on) -- not for deploying on the open internet as-is.
"""

import json
from pathlib import Path

from flask import Flask, render_template, request, redirect, url_for, flash

from .fetch import fetch_structure
from .analyze import analyze_structure
from .domains import annotate_domains

app = Flask(__name__)
app.secret_key = "protein-annotation-toolkit-local-dev"  # fine for local-only use, not for deployment


@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")


@app.route("/structure", methods=["POST"])
def structure_lookup():
    pdb_id = request.form.get("pdb_id", "").strip()
    include_domains = request.form.get("domains") == "on"

    if not pdb_id:
        flash("Enter a PDB ID first.")
        return redirect(url_for("index"))

    return redirect(url_for("structure_view", pdb_id=pdb_id, domains="1" if include_domains else "0"))


@app.route("/structure/<pdb_id>")
def structure_view(pdb_id):
    include_domains = request.args.get("domains") == "1"

    try:
        structure_path = fetch_structure(pdb_id, out_dir="data")
    except Exception as exc:
        flash(f"Couldn't fetch '{pdb_id}': {exc}")
        return redirect(url_for("index"))

    try:
        results = analyze_structure(structure_path)
    except Exception as exc:
        flash(f"Couldn't analyze '{pdb_id}': {exc}")
        return redirect(url_for("index"))

    if include_domains:
        results["domains"] = annotate_domains(pdb_id)

    structure_text = Path(structure_path).read_text()

    return render_template(
        "result.html",
        pdb_id=pdb_id.upper(),
        results=results,
        include_domains=include_domains,
        # json.dumps here, not Jinja's |tojson, so we control ensure_ascii and
        # ship a plain JS string literal for 3Dmol.js's addModel() call.
        structure_json=json.dumps(structure_text),
    )


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Run the protein-annotation-toolkit web viewer")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--debug", action="store_true", default=True)
    args = parser.parse_args()

    print(f"Open http://{args.host}:{args.port} in a browser")
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
