import json
import os
import sys

import pytest

from procsim.__main__ import example_names, main

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_examples_bundled():
    assert {"methanol_loop", "neodymium_bastnaesite", "neodymium_monazite_mtweld",
            "neodymium_eneabba_concentrate"} <= set(example_names())


def test_run_bundled_example_by_name(capsys, tmp_path):
    out = tmp_path / "res.json"
    assert main(["run", "methanol_loop", "-o", str(out)]) == 0
    assert "Recycle" in capsys.readouterr().out
    res = json.loads(out.read_text())
    assert res["ok"] and res["loops"][0]["converged"]


def test_run_json_output(capsys):
    assert main(["run", "copper_concentrator", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["balances"]["mass_in_kg_s"] > 0


def test_run_unknown_reference():
    with pytest.raises(SystemExit):
        main(["run", "no_such_flowsheet"])


def test_examples_copy(tmp_path):
    assert main(["examples", "--copy", str(tmp_path)]) == 0
    assert len(list(tmp_path.glob("*.json"))) == len(example_names())


def test_units_lists_hydrometallurgy(capsys):
    main(["units"])
    out = capsys.readouterr().out
    assert "SXCascade" in out and "Electrolysis" in out


def test_editor_writes_page(tmp_path):
    out = tmp_path / "editor.html"
    assert main(["editor", "--out", str(out), "--no-open"]) == 0
    html = out.read_text(encoding="utf-8")
    assert "<title>Procsim Flowsheet</title>" in html and "ProcSim" in html


def test_routes_json(capsys):
    assert main(["routes", "neodymium_eneabba_concentrate", "--json"]) == 0
    row = json.loads(capsys.readouterr().out)["neodymium_eneabba_concentrate"]
    assert 88 < row["concentrate_to_oxide_pct"] < 93  # Iluka quotes ~90 % refinery recovery


def test_bundled_editor_is_current():
    """procsim/data/procsim.html must be rebuilt (python web/build.py) after editing the web sources or examples."""
    sys.path.insert(0, os.path.join(ROOT, "web"))
    from build import bundle
    with open(os.path.join(ROOT, "procsim", "data", "procsim.html"), encoding="utf-8") as f:
        assert f.read() == bundle(), "run: python web/build.py"
