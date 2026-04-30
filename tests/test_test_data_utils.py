from __future__ import annotations

from pathlib import Path

from gxrender.utils.test_data import (
    find_default_model_file,
    find_ebtel_file,
    find_model_file,
    find_model_files,
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

    (model_dir / "test.chr.h5").write_bytes(b"h5")
    (response_dir / "resp_aia_20251126T153431.sav").write_bytes(b"response")
    (ebtel_dir / "ebtel.sav").write_bytes(b"ebtel")

    monkeypatch.setenv("GXRENDER_TEST_DATA_ROOT", str(raw_root))

    assert find_model_file("test.chr.h5") == model_dir / "test.chr.h5"
    assert find_model_files(".h5") == [model_dir / "test.chr.h5"]
    assert find_default_model_file(".h5") == model_dir / "test.chr.h5"
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
