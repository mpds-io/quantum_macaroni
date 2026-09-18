"""CLI coverage for CRYSTAL bands, DOS, combined inputs, and checkpoint references."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

import main as cli
from quantum_macaroni import calculate_spin_polarized_transport
from quantum_macaroni.core.constants import HTR_TO_EV
from tests.test_crystal_doss import _dat, _fort25
from tests.test_crystal_outp import _ALPHA, _BETA, _HEADER


def _files(tmp_path: Path, jspins: int) -> tuple[Path, Path]:
    """Create paired inputs with intentionally different Fermi energies."""
    outp = tmp_path / "system.outp"
    blocks = "ALPHA ELECTRONS\n" + _ALPHA + "BETA ELECTRONS\n" + _BETA if jspins == 2 else _ALPHA  # noqa: PLR2004
    outp.write_text(_HEADER + blocks)
    doss = tmp_path / "system.DOSS"
    doss.write_text(_dat(jspins))
    return outp, doss


def _run(monkeypatch: pytest.MonkeyPatch, *arguments: str | Path) -> None:
    """Exercise CLI argument parsing, real workflows, and JSON serialization."""
    monkeypatch.setattr(sys, "argv", ["main.py", *(str(argument) for argument in arguments)])
    cli.main()


@pytest.mark.parametrize("jspins", [1, 2])
@pytest.mark.parametrize("fixed_width", [False, True])
def test_dos_cli(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, jspins: int, fixed_width: bool) -> None:
    """Export DOS for either spin layout and either CRYSTAL format."""
    path = tmp_path / "input.data"
    path.write_text(_fort25(jspins - 1) if fixed_width else _dat(jspins))
    output = tmp_path / "dos.json"
    _run(monkeypatch, path, "--output", output)
    result = json.loads(output.read_text())
    assert result["parser"] == "crystal-doss"
    assert result["jspins"] == jspins
    assert np.shape(result["dos"]) == (jspins, 3, 2)
    np.testing.assert_allclose(result["energies"], np.array([-0.1, 0, 0.1]) * HTR_TO_EV, atol=1e-14)
    np.testing.assert_allclose(result["dos"][0], np.array([[1, 10], [2, 20], [3, 30]]) / HTR_TO_EV)
    if jspins == 2:  # noqa: PLR2004
        np.testing.assert_allclose(result["dos"][1], np.array([[4, 40], [5, 50], [6, 60]]) / HTR_TO_EV)
    assert not (tmp_path / "transport_state.npz").exists()


def test_dos_default_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Select DOS export explicitly and write its documented default output."""
    _, doss = _files(tmp_path, 1)
    monkeypatch.chdir(tmp_path)
    _run(monkeypatch, doss, "--parser", "crystal-doss")
    assert json.loads((tmp_path / "dos_results.json").read_text())["jspins"] == 1
    assert not (tmp_path / "transport_state.npz").exists()


@pytest.mark.parametrize("jspins", [1, 2])
@pytest.mark.parametrize("positional_dos", [False, True])
def test_combined_cli(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, jspins: int, positional_dos: bool) -> None:
    """Run band transport and include DOS under the same JSON root."""
    outp, doss = _files(tmp_path, jspins)
    output = tmp_path / "transport.json"
    companion = [str(doss)] if positional_dos else ["--dos-file", str(doss)]
    _run(
        monkeypatch,
        outp,
        *companion,
        "--kmesh",
        "3",
        "3",
        "3",
        "--lr-ratio",
        "2",
        "--no-checkpoint",
        "--output",
        output,
    )
    result = json.loads(output.read_text())
    assert result["meta"]["parser"] == "crystal-outp"
    assert result["meta"]["jspins"] == jspins
    assert result["meta"]["fermi_energy"] == pytest.approx(-0.2 * HTR_TO_EV)
    assert result["meta"]["fermi_source"] == "doss"
    assert result["dos"]["jspins"] == jspins
    assert np.shape(result["dos"]["dos"]) == (jspins, 3, 2)
    assert np.isfinite(result["0.0"]["300.0"]["sigma"]).all()


def test_keep_outp_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Honor an explicit choice to keep the outp Fermi energy."""
    outp, doss = _files(tmp_path, 1)
    output = tmp_path / "transport.json"
    _run(
        monkeypatch,
        outp,
        "--doss",
        doss,
        "--fermi-source",
        "outp",
        "--kmesh",
        "3",
        "3",
        "3",
        "--lr-ratio",
        "2",
        "--no-checkpoint",
        "--output",
        output,
    )
    result = json.loads(output.read_text())
    assert result["meta"]["fermi_energy"] == pytest.approx(-0.15 * HTR_TO_EV)
    assert result["meta"]["fermi_source"] == "outp"
    assert result["dos"]["fermi_energy"] == pytest.approx(-0.2 * HTR_TO_EV)


def test_outp_cli(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Run outp alone with automatic parser detection."""
    outp, _ = _files(tmp_path, 1)
    output = tmp_path / "transport.json"
    _run(monkeypatch, outp, "--kmesh", "3", "3", "3", "--lr-ratio", "2", "--no-checkpoint", "--output", output)
    result = json.loads(output.read_text())
    assert result["meta"]["parser"] == "crystal-outp"
    assert result["meta"]["fermi_energy"] == pytest.approx(-0.15 * HTR_TO_EV)
    assert "dos" not in result


def test_dos_reference_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Invalidate transport when DOS Fermi changes, then reuse unchanged transport."""
    outp, doss = _files(tmp_path, 1)
    output = tmp_path / "transport.json"
    checkpoint = tmp_path / "transport.npz"
    arguments = [
        outp,
        doss,
        "--kmesh",
        "3",
        "3",
        "3",
        "--lr-ratio",
        "2",
        "--checkpoint",
        checkpoint,
        "--output",
        output,
    ]
    _run(monkeypatch, *arguments)
    capsys.readouterr()
    doss.write_text(_dat(1).replace("-0.20000D+00", "-0.18000D+00"))
    _run(monkeypatch, *arguments)
    result = json.loads(output.read_text())
    assert result["meta"]["fermi_energy"] == pytest.approx(-0.18 * HTR_TO_EV)
    assert "Loaded completed transport result" not in capsys.readouterr().out
    doss.write_text(doss.read_text().replace("1.0000D+00", "1.0000D+02"))
    _run(monkeypatch, *arguments)
    result = json.loads(output.read_text())
    assert result["dos"]["dos"][0][0][0] == pytest.approx(100 / HTR_TO_EV)
    assert "Loaded completed transport result" in capsys.readouterr().out


@pytest.mark.parametrize("invalid", ["spin", "dos-primary", "duplicate", "missing", "bad-dos", "xml"])
def test_cli_input_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, invalid: str) -> None:
    """Reject invalid input pairs before writing output or checkpoints."""
    outp, doss = _files(tmp_path, 1)
    output = tmp_path / "transport.json"
    arguments: list[str | Path] = [outp, "--dos-file", doss]
    if invalid == "spin":
        doss.write_text(_dat(2))
    elif invalid == "dos-primary":
        arguments[0] = doss
    elif invalid == "duplicate":
        arguments.insert(1, doss)
    elif invalid == "missing":
        arguments[0] = tmp_path / "missing.outp"
    elif invalid == "bad-dos":
        doss.write_text("bad DOS data")
    else:
        arguments.extend(["--parser", "fleur-outxml"])
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as exc:
        _run(monkeypatch, *arguments, "--output", output)
    assert exc.value.code == 2  # noqa: PLR2004
    assert not output.exists()
    assert not (tmp_path / "transport_state.npz").exists()


def test_xml_detection() -> None:
    """Preserve FLEUR XML support under the automatic selector."""
    path = Path(__file__).resolve().parents[1] / "examples" / "PbTe-nospin" / "out-nospin.xml"
    assert cli._input_parser(str(path)) == "fleur-outxml"


@pytest.mark.parametrize("fermi_energy", [float("nan"), float("inf")])
def test_invalid_reference(tmp_path: Path, fermi_energy: float) -> None:
    """Validate overrides before parsing input or creating a checkpoint."""
    with pytest.raises(ValueError, match="fermi_energy must be finite"):
        calculate_spin_polarized_transport(tmp_path / "missing.outp", fermi_energy=fermi_energy)
