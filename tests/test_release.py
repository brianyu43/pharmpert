from pathlib import Path

import pandas as pd
import pytest

from yakseopdong.release import sha256
from yakseopdong.release_validation import _validate_manifest


def test_sha256_reads_file_in_binary_chunks(tmp_path: Path) -> None:
    path = tmp_path / "artifact.bin"
    path.write_bytes(b"pharmpert-release\n")
    assert sha256(path) == "1a2052b3b9132d94bd70c0494fdc2a9a9c714da8e62bd91c2c88354599c342db"


def test_manifest_validation_checks_hash_rows_and_bytes(tmp_path: Path) -> None:
    path = tmp_path / "table.csv"
    path.write_text("value\n1\n2\n", encoding="utf-8")
    manifest = pd.DataFrame(
        [
            {
                "path": "table.csv",
                "rows": 2,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        ]
    )
    _validate_manifest(tmp_path, manifest, is_table=True)
    manifest.loc[0, "rows"] = 3
    with pytest.raises(ValueError, match="row count"):
        _validate_manifest(tmp_path, manifest, is_table=True)
