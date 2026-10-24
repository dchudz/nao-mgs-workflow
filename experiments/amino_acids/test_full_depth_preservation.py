"""Private evidence must be verified before any local copy is evicted."""
import io,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch

class PreservationTest(unittest.TestCase):
    def test_eviction_requires_identical_private_copy(self):
        source=(Path(__file__).with_name('verify_full_depth_backup.py')).read_text()
        prefix='s3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/'
        for corrupt in [True,False]:
            with self.subTest(corrupt=corrupt),tempfile.TemporaryDirectory() as temporary:
                root=Path(temporary);sample='SRRtest';d=root/'e2e'/sample
                plan=root/'experiments/amino_acids/end-to-end-full-depth';plan.mkdir(parents=True)
                (plan/'sample-selection.json').write_text(json.dumps({'initial_samples':[{'run_accession':sample}]}))
                d.mkdir(parents=True);(d/'complete.json').write_text(json.dumps({'all_four_outputs_completed':True,'source_commit':'a'*40}))
                for stage in ['PROTEIN_ALIGN','PROTEIN_CANDIDATES']:
                    target=d/'protein-run/published'/stage/'protein-hits.tsv';target.parent.mkdir(parents=True);target.write_bytes(b'query\treference\n')
                work=d/'protein-run/work/aa/bb';work.mkdir(parents=True);(work/'metrics.json').write_text('{}')
                (d/'protein-run/task-diagnostics.tar.gz').write_bytes(b'diagnostics')
                final=d/'protein_variant-downstream/output/results_downstream/validation_hits.tsv.gz';final.parent.mkdir(parents=True);final.write_bytes(b'final scientific result')
                raw=root/'data'/sample;raw.mkdir(parents=True)
                for name in [sample,'R1.fastq.gz','R2.fastq.gz']:(raw/name).write_bytes(b'public raw input')
                remote={prefix+'samples/'+sample+'/'+str(f.relative_to(d)):f.read_bytes() for f in d.rglob('*') if f.is_file() and 'work' not in f.relative_to(d).parts}
                def run(cmd,**kwargs):
                    self.assertEqual(cmd[5:7],['s3','cp']);remote[cmd[8]]=Path(cmd[7]).read_bytes()
                class Download:
                    def __init__(self,cmd,**kwargs):
                        data=remote[cmd[7]]
                        if corrupt and cmd[7].endswith('protein-hits.tsv'):data+=b'corrupted'
                        self.stdout=io.BytesIO(data)
                    def wait(self,timeout=None):return 0
                    def poll(self):return 0
                program=source.replace("R=Path('/home/ec2-user/experiment')",f'R=Path({str(root)!r})')
                argv=['verify',sample,'--release-scratch','--release-published-hsps']
                with patch.object(sys,'argv',argv),patch('subprocess.check_output',return_value=json.dumps({'Account':'839993661825'})),patch('subprocess.run',side_effect=run),patch('subprocess.Popen',Download),patch('builtins.print'):
                    if corrupt:
                        with self.assertRaisesRegex(AssertionError,'Independent readback differs'):exec(compile(program,'verify','exec'),{'__name__':'__main__'})
                        self.assertTrue(work.exists());self.assertTrue((raw/'R1.fastq.gz').exists())
                        self.assertEqual(len(list(d.glob('protein-run/published/*/protein-hits.tsv'))),2)
                    else:
                        exec(compile(program,'verify','exec'),{'__name__':'__main__'})
                        self.assertFalse(work.exists());self.assertFalse((raw/'R1.fastq.gz').exists());self.assertTrue(final.exists())
                        receipt=json.loads((d/'backup-verification.json').read_text());manifest=json.loads((d/'file-manifest.json').read_text())
                        self.assertEqual(len(receipt['locally_evicted_hsps']),2)
                        for item in receipt['locally_evicted_hsps']:
                            self.assertFalse((d/item['path']).exists());self.assertEqual(remote[item['s3_uri']],b'query\treference\n')
                            self.assertTrue(any(m['path']==item['path'] and m['sha256']==item['sha256'] for m in manifest))
                        with self.assertRaisesRegex(AssertionError,'Restore the recorded S3 HSPs'):exec(compile(program,'verify','exec'),{'__name__':'__main__'})

if __name__=='__main__':unittest.main()
