# Execution halt, preservation and spending

Execution was halted by the assistant, not requested by David. The previously accepted expansion crossed the assistant's boundary on improving discovery of divergent human-infecting viruses. That acceptance was mistaken. The halt was not triggered by cost, experimental output, or a correctness failure. No further pipeline computation is scheduled.

Seventeen paired comparisons are complete, accepted, individually committed, and protected by full-readback-verified private S3 archives. Forty-three planned comparisons remain unfinished. No unfinished unit is reported as zero. The frozen known-origin artificial control passed before expansion; its separate report remains the authoritative control record.

## Worker and data disposition

| Worker | Verified disposition | Retained volumes |
| --- | --- | --- |
| primary (i-02d118993a6f4989a) | stopped at 2026-10-24T13:14:32.266264+00:00 | vol-0272ade891333af83 |
| spot2 (i-0ee9da3f09d197be6) | terminated at 2026-10-24T13:14:32.042477+00:00 | vol-05ca1445d265f5325, vol-03f6a458f0812c7b4 |

The local foreground supervisor was terminated after its worker processes were suspended. Worker logs and completion inventories were uploaded and fully read back before shutdown. The primary instance remains stopped. The one-time Spot instance was terminated with DeleteOnTermination=false on both attached disks; both disks remain available. Incomplete filesystem output is retained, but no claim is made that an interrupted task is complete or restartable without verification.

The earlier interrupted Spot worker's root vol-004c4ffab08c46dbd and scratch vol-06f002b044d0aeb35 also remain retained. The original reference snapshot snap-0c13ccf01c762eb8a is completed and retained. No volume, snapshot, reference, or result was deleted.

## Private preservation locations

- primary: `s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/preservation-halt/primary-halt-metadata.tar.gz`, 6,657,413 bytes; SHA-256 `25e2199519e89c1f3b1ff05df33dbd2d914b0289832e6071346c70505994b037`; full readback verified.
- spot2: `s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/preservation-halt/spot2-halt-metadata.tar.gz`, 6,774,530 bytes; SHA-256 `2aef2db0e0ee6785eb249559a332f9fccf43a7b526c65de934c0a1de95111d33`; full readback verified.

Completed full evidence archives remain under `s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/comparisons/<label>/completed-evidence.tar.gz`. Each accepted result directory contains its exact location and checksum. Selected-input archives have separate readback receipts. The two halt archives contain metadata and logs; they are not substitutes for complete evidence archives or retained disks.

## Cost record

At 2026-10-24T13:14:54.515287+00:00, the conservative new-compute ceiling was **$47.94**, and modeled study-period EBS accrual was **$7.26**. All compute intervals are now closed using verified terminal-state timestamps. These are modeled bounds, not invoices.

The five retained volumes cost about **$753.44 per 30-day month ($25.11/day)** under the existing rate model. Storage continues after compute stops. S3 currently holds 92,788,830,328 bytes across 3,539 objects for this study prefix. Snapshot incremental storage, S3 age/requests, transfer and taxes remain unpriced in this model. No all-in billed total is claimed. The $500 budget and its reserves were not exhausted.

The complete delayed account-wide Cost Explorer response is saved in spend-halted.json. It includes overlapping work and must not be added to the compute/storage models as a separate study expense. Measured per-stage wall times remain in results/aggregate/compute-split.md; acquisition, preservation and idle intervals are distinct. The CloudWatch snapshot includes coarse five-minute CPU and disk metrics, not stage-level CPU or memory measurements.

## Interpretation

The completed comparisons measure output differences under the frozen pipeline policies. They do not establish sensitivity, precision, validated new species calls, or a preferred deployment method. Downstream settings match, but upstream scoring thresholds and false-positive operating points do not. The incomplete cohort cannot support the originally planned full-cohort comparison. No candidate-specific reconstruction or new search was added during preservation.
