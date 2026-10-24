"""Reject comparison acceptance when the known-answer control is stale or absent."""
import hashlib
import json
import subprocess
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_control(source, diamond, image, control_dir=None):
    source=Path(source)
    control=Path(control_dir) if control_dir else source/'experiments/paired_pipeline_comparison'
    fixture=control/'synthetic-control-fixture'
    receipt=control/'synthetic-control-results/attempt-3/result.json'
    result=json.loads(receipt.read_text())
    if result.get('passed') is not True:
        raise ValueError('Synthetic control did not pass')
    checks={
        fixture/'expected.json':result['expectation_sha256'],
        fixture/'manifest.json':result['fixture_manifest_sha256'],
        control/'check_synthetic_adapter.py':result['runner_sha256'],
        **{source/name:value for name,value in result['source_hashes'].items()},
    }
    manifest=json.loads((fixture/'manifest.json').read_text())
    checks.update({fixture/name:value for name,value in manifest.items()})
    for path,expected in checks.items():
        if digest(path)!=expected:
            raise ValueError(f'Synthetic control is stale for {path}')
    expected=json.loads((fixture/'expected.json').read_text())
    for field in ('positive_pairs','exact_unique_pairs','alignment_exemplars','joint_exemplars'):
        if result[field]!=expected[field]:
            raise ValueError(f'Synthetic control count differs: {field}')
    if result['geometry_cases']!=36 or result['negative_pairs_excluded']!=len(expected['negative_ids']):
        raise ValueError('Synthetic control coverage is incomplete')
    if set(result['rejected_error_injections'])!={'prim_align_ref_start','aligner_taxid_lca','prim_align_dup_exemplar'}:
        raise ValueError('Synthetic control did not reject the injected errors')
    diamond_hash=next(value for name,value in result['binary_hashes'].items() if Path(name).name=='diamond')
    if digest(diamond)!=diamond_hash:
        raise ValueError('DIAMOND differs from the tested binary')
    image_id=subprocess.check_output(['docker','image','inspect',image,'--format','{{.Id}}'],text=True).strip()
    if image_id!=result['duplicate_container_id']:
        raise ValueError('Duplicate container differs from the tested image')
    return {'passed':True,'receipt_sha256':digest(receipt),'duplicate_container_id':image_id,
            'scope':'Known-origin synthetic adapter and duplicate handling only; other comparison integrity checks still apply.'}
