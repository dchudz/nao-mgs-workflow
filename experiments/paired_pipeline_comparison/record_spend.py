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
VOLUMES = ['vol-0272ade891333af83', 'vol-004c4ffab08c46dbd', 'vol-05ca1445d265f5325', 'vol-06f002b044d0aeb35', 'vol-03f6a458f0812c7b4']
AWS = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']

def get(*args):
    return json.loads(subprocess.check_output(AWS + list(args), text=True))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--intervals', type=Path, required=True)
    a = p.parse_args()
    assert get('sts', 'get-caller-identity')['Account'] == '839993661825'
    intervals = json.loads(a.intervals.read_text())['intervals']
    active_ids = sorted({x['instance_id'] for x in intervals if x['ended_utc_upper'] is None})
    assert set(active_ids) <= set(WORKERS)
    # EC2 eventually omits terminated instances. Closed, evidenced intervals
    # remain billable history even after their live instance record disappears.
    instances = get('ec2', 'describe-instances', '--instance-ids', *active_ids) if active_ids else {'Reservations': []}
    now = dt.datetime.now(dt.timezone.utc)
    inventory = {x['InstanceId']:x for r in instances['Reservations'] for x in r['Instances']}
    rows = []
    for instance_id, (worker, original_start, rate) in WORKERS.items():
        periods = sorted((x for x in intervals if x['instance_id']==instance_id), key=lambda x:x['started_utc'])
        assert periods and periods[0]['started_utc']==original_start
        opened = [x for x in periods if x['ended_utc_upper'] is None]
        state = inventory.get(instance_id, {}).get('State', {}).get('Name', 'closed_interval_no_live_query')
        assert len(opened)==(1 if state in {'running','pending','stopping','shutting-down'} else 0)
        if opened:
            assert dt.datetime.fromisoformat(inventory[instance_id]['LaunchTime'])==dt.datetime.fromisoformat(opened[0]['started_utc'])
        last = None
        for period in periods:
            start = dt.datetime.fromisoformat(period['started_utc'])
            end = dt.datetime.fromisoformat(period['ended_utc_upper']) if period['ended_utc_upper'] else now
            assert start<=end<=now and (last is None or start>=last)
            last=end;hours=(end-start).total_seconds()/3600
            rows.append(dict(worker=worker,instance_id=instance_id,started_utc=start.isoformat(),
                             ended_utc_upper=period['ended_utc_upper'],observed_state=state,
                             observed_utc=now.isoformat(),hours=hours,hourly_compute_ceiling_USD=rate,
                             compute_ceiling_USD=hours*rate,evidence=period.get('evidence')))
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
    billing = get('ce', 'get-cost-and-usage', '--time-period', 'Start=2026-10-23,End='+(now.date()+dt.timedelta(days=1)).isoformat(),
                  '--granularity', 'DAILY', '--metrics', 'UnblendedCost', '--group-by', 'Type=DIMENSION,Key=SERVICE')
    result = dict(observed_utc=now.isoformat(), account='839993661825', workers=rows, volumes=storage,
                  elapsed_compute_ceiling_USD=sum(r['compute_ceiling_USD'] for r in rows),
                  new_compute_allowance_USD=200, study_EBS_accrual_model_USD=sum(r['accrued_model_USD'] for r in storage),
                  study_S3_current_objects=len(objects), study_S3_current_bytes=sum(o['Size'] for o in objects),
                  cost_explorer_account_wide=billing,
                  unpriced=['Snapshot incremental storage', 'S3 storage age and requests', 'Data transfer', 'Taxes'],
                  limits=['Compute is a sum of disjoint observed running intervals at rate ceilings, not an invoice. Closed intervals use a conservative verified-stop time.',
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
