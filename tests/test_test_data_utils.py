from __future__ import annotations

from pathlib import Path

from gxrender.utils.test_data import (
    find_default_model_file,
    find_ebtel_file,
    find_model_file,
    find_model_files,
    find_model_loader_parity_files,
    find_response_file,
)


def test_test_data_locator_uses_external_raw_layout(monkeypatch, tmp_path: Path) -> None:
    raw_root = tmp_path / "raw"
    model_dir = raw_root / "models" / "models_20251126T153431"
    response_dir = raw_root / "responses" / "20251126T153431"
    ebtel_dir = raw_root / "ebtel" / "ebtel_gxsimulator_euv"

    model_dir.mkdir(parents=True)
    response_dir.mkdir(parents=True)
    ebtel_dir.mkdir(parents=True)

    model_file = model_dir / "hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.GEN.CHR.h5"
    model_file.write_bytes(b"h5")
    (response_dir / "resp_aia_20251126T153431.sav").write_bytes(b"response")
    (ebtel_dir / "ebtel.sav").write_bytes(b"ebtel")

    monkeypatch.setenv("GXRENDER_TEST_DATA_ROOT", str(raw_root))

    assert find_model_file(model_file.name) == model_file
    assert find_model_files(".h5") == [model_file]
    assert find_default_model_file(".h5") == model_file
    assert find_response_file("aia") == response_dir / "resp_aia_20251126T153431.sav"
    assert find_ebtel_file("ebtel.sav") == ebtel_dir / "ebtel.sav"


def test_default_model_locator_falls_back_to_later_existing_root(monkeypatch, tmp_path: Path) -> None:
    empty_raw = tmp_path / "empty" / "raw"
    populated_raw = tmp_path / "populated" / "raw"
    model_dir = populated_raw / "models" / "models_20201126T195831"
    empty_raw.mkdir(parents=True)
    model_dir.mkdir(parents=True)
    model_path = model_dir / "renamed.chr.h5"
    model_path.write_bytes(b"h5")

    monkeypatch.setenv("GXRENDER_TEST_DATA_ROOT", str(empty_raw))
    monkeypatch.setenv("GXIMAGECOMPUTING_TEST_DATA_ROOT", str(populated_raw))

    assert find_model_files(".h5") == [model_path]
    assert find_default_model_file(".h5") == model_path


def test_model_locator_prefers_default_h5_model_over_loader_parity_clone(monkeypatch, tmp_path: Path) -> None:
    raw_root = tmp_path / "raw"
    model_dir = raw_root / "models" / "models_20201126T195831"
    parity_dir = raw_root / "models" / "model_loader_parity_20201126T195831"
    model_dir.mkdir(parents=True)
    parity_dir.mkdir(parents=True)

    default_h5 = model_dir / "hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.GEN.CHR.h5"
    sav_path = parity_dir / "hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.CHR.sav"
    clone_h5 = parity_dir / "hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.CHR.clone.h5"
    default_h5.write_bytes(b"h5")
    sav_path.write_bytes(b"sav")
    clone_h5.write_bytes(b"clone")

    monkeypatch.setenv("GXRENDER_TEST_DATA_ROOT", str(raw_root))

    assert find_default_model_file(".h5") == default_h5
    assert find_model_loader_parity_files() == (sav_path, clone_h5)
