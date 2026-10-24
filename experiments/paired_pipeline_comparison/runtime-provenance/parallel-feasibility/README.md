# Parallel feasibility execution

The first paired comparison completed, including its assignment-output integrity
checks. Its 52,297,836-byte evidence archive was uploaded to private S3 and fully
read back; SHA-256 was identical. Both arms have an explicitly unavailable shared
Bracken profile, so this does not establish a successful abundance comparison.

The remaining nine feasibility runs are assigned before their results are
inspected: three on the retained primary worker, three on each of two private
Spot replicas. The two prespecified nested 25-million-pair cases are on separate
replicas. All three use source c77cdde08b03834a9c1f54f3d80b0c9668bf9a83.

The replicas use the retained private snapshot and the same verified sampler,
DIAMOND binary, duplicate container and source files. SSH host keys were checked
against EC2 console output. Both replicas passed the synthetic-control identity
gate before their batches started. This verifies identity with the tested
implementation; it is not a new accuracy benchmark.

Both Spot requests are one-time requests capped at $2/hour. Their 96 vCPUs fit
the existing 128-vCPU Spot quota. The primary uses the existing 48-vCPU on-demand
quota. Shutdown guards end by 2026-10-25 05:07 UTC, with root volumes retained.
The primary stops; the Spot instances terminate. The shared maximum allocation
is $192 of the $200 new-compute allowance, within the existing $500 overall cap.
This is a spending safeguard, not a user-imposed deadline. Completed workers
should end sooner after backups are verified.

At the recorded observation, elapsed compute was bounded by $6.06. This is not
an invoice and excludes storage and transfer. Cost Explorer still had no posted
October 24 charges. The existing budget keeps $150 reserved for billing lag and
noncompute expenses; persistent replica volumes need cleanup after preservation.

No calling rules, alignment parameters, reference content, or validation
settings changed for parallel execution. Final reports must distinguish output
differences from accuracy, and shared-QC component times from independent
standalone pipeline timings.
