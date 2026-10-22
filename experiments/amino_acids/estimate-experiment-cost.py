"""Model owned-worker lifetime cost; never substitute this model for delayed billing."""
import argparse, datetime, json
from pathlib import Path

r = Path(__file__).resolve().parent
p = argparse.ArgumentParser()
p.add_argument('--as-of', default=datetime.datetime.now(datetime.timezone.utc).isoformat())
p.add_argument('--retirement-receipts', help='JSON mapping instance IDs to confirmed termination UTC timestamps')
a = p.parse_args()
retired = json.loads(Path(a.retirement_receipts).read_text()) if a.retirement_receipts else {}
now = datetime.datetime.fromisoformat(a.as_of.replace('Z', '+00:00'))
def dt(x):
    return datetime.datetime.fromisoformat(x.replace('Z', '+00:00'))
def hours(start, end):
    return max(0, (min(dt(end), now) - dt(start)).total_seconds() / 3600)
primary = 'i-061f366bd55def34c'
nr = 'i-0103f018916b6c9e5'
screen = 'i-0dd39a46c4efb6c4c'
screen_end = '2026-10-22T05:22:19.224178Z'
intervals = [
    (primary, 'x2iedn.4xlarge', '2026-10-21T15:24:12Z', retired.get(primary, a.as_of), 3.3345),
    (screen, 'c7i.4xlarge', '2026-10-21T16:19:35Z', '2026-10-21T22:00:13Z', .714),
    (screen, 'c7i.2xlarge', '2026-10-21T22:00:13Z', screen_end, .357),
    (nr, 'r7i.4xlarge', '2026-10-21T17:11:13Z', '2026-10-21T22:00:13Z', 1.0584),
    (nr, 'm8azn.6xlarge', '2026-10-21T22:00:13Z', '2026-10-22T04:06:26Z', 2.47752),
    (nr, 'r7a.8xlarge', '2026-10-22T05:22:36.767549Z', retired.get(nr, a.as_of), 2.4344),
    ('i-048be724ac0462256', 'm7i.16xlarge Spot', '2026-10-21T18:06:37Z', '2026-10-21T18:50:53Z', 1.5461),
    ('i-0f577952bbbd27d36', 'c6a.24xlarge Spot', '2026-10-21T21:32:32Z', '2026-10-21T21:51:41Z', 1.1061),
]
compute = [{'instance': i, 'type': t, 'start': s, 'end': e, 'modeled_hours': hours(s, e), 'USD_per_hour': rate, 'modeled_USD': hours(s, e) * rate} for i, t, s, e, rate in intervals]
volumes = [
    ('vol-098cdf9d650bd4a0d', primary, '2026-10-21T15:24:13.296Z', retired.get(primary, a.as_of), 2048, 6000, 1000),
    ('vol-07d4e85519abc9f76', nr, '2026-10-21T17:11:14.288Z', retired.get(nr, a.as_of), 1280, 6000, 1000),
    ('vol-0461a97dae04945e3', screen, '2026-10-21T16:19:36.008Z', screen_end, 256, 3000, 500),
    ('vol-0aad5b7d9f190c48d', 'i-048be724ac0462256', '2026-10-21T18:06:37Z', '2026-10-21T18:50:53Z', 512, 6000, 1000),
    ('vol-0533d3e9b96594f40', 'i-0f577952bbbd27d36', '2026-10-21T21:32:32Z', '2026-10-21T21:51:41Z', 1280, 6000, 1000),
]
storage = []
for volume, instance, start, end, gib, iops, throughput in volumes:
    rate = (gib * .08 + max(0, iops - 3000) * .005 + max(0, throughput - 125) * .04) / 720
    storage.append({'volume': volume, 'instance': instance, 'start': start, 'end': end, 'GiB': gib, 'IOPS': iops, 'MiB_per_second': throughput, 'maximum_configuration_modeled_USD': hours(start, end) * rate})
result = {
    'as_of_utc': a.as_of,
    'account': '839993661825',
    'compute_intervals': compute,
    'volume_intervals': storage,
    'modeled_EC2_compute_USD': sum(x['modeled_USD'] for x in compute),
    'modeled_EBS_maximum_configuration_USD': sum(x['maximum_configuration_modeled_USD'] for x in storage),
    'modeled_public_IPv4_USD': sum(x['modeled_hours'] for x in compute) * .005,
    'interpretation': [
        'Lifetime allocation includes reference staging, failed and interrupted attempts, concurrent searches, analyses, and idle time. Individual search allocations must not be added to this ledger.',
        'Short resize stop gaps before 22:00 are conservatively counted as running; the longer NR stop from 04:06:26 to 05:22:36 is excluded from compute but included in EBS.',
        'EBS uses the maximum observed configuration for the entire lifetime and a conservative 720-hour month. The primary volume initially had lower IOPS and throughput.',
        'Historical Spot rates are verified in primary-closure/spot-price-*.json. Runtime times rate does not account for potential AWS interruption credits.',
        'S3 storage, requests, data transfer, taxes, credits and unrelated account usage are excluded. Cost Explorer and retained-storage estimates must be shown separately.',
        'This is a resource-based estimate, not an invoice or a claim that Cost Explorer is current.',
    ],
}
result['modeled_owned_compute_EBS_IPv4_total_USD'] = sum(result[k] for k in ['modeled_EC2_compute_USD', 'modeled_EBS_maximum_configuration_USD', 'modeled_public_IPv4_USD'])
out = r / 'results/cost-ledger.json'
out.write_text(json.dumps(result, indent=2))
print(json.dumps({k: v for k, v in result.items() if k not in ['compute_intervals', 'volume_intervals', 'interpretation']}))
