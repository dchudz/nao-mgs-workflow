"""Record observed resource time separately from delayed account-wide billing."""
import argparse
import datetime as dt
import json
import subprocess
from pathlib import Path

WORKERS = {
    'i-02d118993a6f4989a': ('primary', '2026-10-24T05:07:00+00:00', 4),
    'i-0f5065e3b4f67cfbc': ('spot1', '2026-10-24T06:10:29+00:00', 2),
    'i-0ee9da3f09d197be6': ('spot2', '2026-10-24T06:11:40+00:00', 2),
}
VOLUMES = ['vol-0272ade891333af83', 'vol-004c4ffab08c46dbd', 'vol-05ca1445d265f5325']
AWS = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']

def get(*args):
    return json.loads(subprocess.check_output(AWS + list(args), text=True))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    assert get('sts', 'get-caller-identity')['Account'] == '839993661825'
    instances = get('ec2', 'describe-instances', '--instance-ids', *WORKERS)
    now = dt.datetime.now(dt.timezone.utc)
    rows = []
    for reservation in instances['Reservations']:
        for instance in reservation['Instances']:
            name, start, rate = WORKERS[instance['InstanceId']]
            # A restart needs a new interval, not an overwritten start time.
            assert dt.datetime.fromisoformat(instance['LaunchTime']) == dt.datetime.fromisoformat(start)
            assert instance['State']['Name'] == 'running', 'Record the stopped interval before refreshing this ledger'
            hours = (now - dt.datetime.fromisoformat(start)).total_seconds() / 3600
            rows.append(dict(worker=name, instance_id=instance['InstanceId'], started_utc=start,
                             state='running', observed_utc=now.isoformat(), hours=hours,
                             hourly_compute_ceiling_USD=rate, compute_ceiling_USD=hours*rate))
    volumes = get('ec2', 'describe-volumes', '--volume-ids', *VOLUMES)['Volumes']
    storage = []
    study_start = dt.datetime.fromisoformat(WORKERS['i-02d118993a6f4989a'][1])
    for volume in volumes:
        assert volume['VolumeType'] == 'gp3'
        start = max(study_start, dt.datetime.fromisoformat(volume['CreateTime']))
        hours = (now-start).total_seconds()/3600
        # Previously sourced US-East model, not an invoice. 30-day month.
        monthly = volume['Size']*.08 + max(0, volume['Iops']-3000)*.005 + max(0, volume['Throughput']-125)*.04
        storage.append(dict(volume_id=volume['VolumeId'], size_GiB=volume['Size'], iops=volume['Iops'],
                            throughput_MiB_per_s=volume['Throughput'], study_hours=hours,
                            modeled_USD_per_30day_month=monthly, accrued_model_USD=monthly*hours/720))
    objects = get('s3api', 'list-objects-v2', '--bucket', 'pipeline-lab-839993661825-us-east-1',
                  '--prefix', 'runs/paired-comparison-20261024/').get('Contents', [])
    billing = get('ce', 'get-cost-and-usage', '--time-period', 'Start=2026-10-23,End=2026-10-25',
                  '--granularity', 'DAILY', '--metrics', 'UnblendedCost', '--group-by', 'Type=DIMENSION,Key=SERVICE')
    result = dict(observed_utc=now.isoformat(), account='839993661825', workers=rows, volumes=storage,
                  elapsed_compute_ceiling_USD=sum(r['compute_ceiling_USD'] for r in rows),
                  new_compute_allowance_USD=200, study_EBS_accrual_model_USD=sum(r['accrued_model_USD'] for r in storage),
                  study_S3_current_objects=len(objects), study_S3_current_bytes=sum(o['Size'] for o in objects),
                  cost_explorer_account_wide=billing,
                  unpriced=['Snapshot incremental storage', 'S3 storage age and requests', 'Data transfer', 'Taxes'],
                  limits=['Compute is an elapsed-time ceiling, not an invoice.',
                          'EBS model excludes pre-study retention and uses a 30-day month.',
                          'S3 bytes are a current inventory, not byte-hours or a charge.',
                          'Cost Explorer is account-wide and delayed; do not add it to these overlapping models.',
                          'No posted charge for a day means unposted, not free.',
                          'The existing 150 USD billing-lag/noncompute reserve remains separate.'])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('observed_utc','elapsed_compute_ceiling_USD','study_EBS_accrual_model_USD','study_S3_current_objects','study_S3_current_bytes')}))

if __name__ == '__main__':
    main()
