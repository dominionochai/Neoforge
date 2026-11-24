from neoforge.evidence import FeatureValue, EvidenceKind, measurement, missing, prediction


def test_feature_value_serialisation():
    fv = prediction(25.8, "mhcflurry-2.4.0")
    d = fv.as_dict()
    assert d["kind"] == "predicted" and d["value"] == 25.8 and d["source"].startswith("mhcflurry")


def test_missing_never_carries_value():
    fv = missing("truth VCF has no tumour VAF field")
    assert fv.kind is EvidenceKind.MISSING and fv.value is None and fv.note
    assert measurement(1.2, "ccle").kind is EvidenceKind.MEASURED
