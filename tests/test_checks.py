from pathlib import Path

from yakseopdong.checks import smoke_report


def test_stage0_smoke_report_is_ok() -> None:
    root = Path(__file__).resolve().parents[1]
    report = smoke_report(root)
    assert report["ok"] is True
    assert report["missing_paths"] == []
    assert report["project"]["generalization_unit"] == "cell_line"


def test_blueprint_copy_is_in_sync() -> None:
    root = Path(__file__).resolve().parents[1]
    assert (root / "yak_seopdong_project_blueprint.md").read_bytes() == (
        root / "docs/project_blueprint.md"
    ).read_bytes()
