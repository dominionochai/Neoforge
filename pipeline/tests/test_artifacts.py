import hashlib
import json

from neoforge.artifacts import ArtifactStore, StageRun


def test_artifact_write_and_load_roundtrip(tmp_path):
    store = ArtifactStore(tmp_path / "work")
    run = StageRun(stage="ingest", case="demo", inputs={"a.vcf.gz": "deadbeef"})
    run.record_software({"mhcflurry": "2.4.0"})
    artifact_id, digest = store.write(run, "ingest", {"n": 42})

    stored = json.loads((tmp_path / "work" / "demo" / f"ingest.{artifact_id}.json").read_bytes())
    assert stored == {"n": 42}

    manifest = store.read_manifest("demo")
    assert manifest["stages"]["ingest"]["artifact_id"] == artifact_id
    assert manifest["stages"]["ingest"]["sha256"] == hashlib.sha256(
        json.dumps({"n": 42}, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert store.load("demo", "ingest") == {"n": 42}


def test_artifact_bytes_are_hash_stable(tmp_path):
    store = ArtifactStore(tmp_path)
    run = StageRun(stage="variants", case="demo")
    a_id, a_sha = store.write(run, "variants", {"b": 1, "a": 2})
    b_id, b_sha = store.write(StageRun(stage="variants", case="demo"), "variants", {"a": 2, "b": 1})
    assert a_id == b_id and a_sha == b_sha


def test_load_missing_artifact_raises(tmp_path):
    store = ArtifactStore(tmp_path)
    import pytest
    with pytest.raises(FileNotFoundError):
        store.load("case", "nonexistent")
