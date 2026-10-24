# Resume of the frozen paired comparison

David explicitly instructed resumption of the remaining 43 comparisons under the
existing $500 budget and quotas. The earlier halt was an incorrect scope decision;
it was not caused by the data, correctness control, or budget. The historical halt
record remains preserved and is superseded by execution-status.json for live status.

The 17 completed comparisons are excluded. The primary worker receives 21 remaining
units across 19 runs, and the replacement Spot worker receives 22 across 19 runs.
Only orchestration changes: explicit remaining depth lists prevent rerunning an
already-complete shallow comparison when its deeper comparison remains unfinished.
The scientific driver, sampling, references, searches, duplicate rules and downstream
settings remain pinned to c77cdde08b03834a9c1f54f3d80b0c9668bf9a83. Both workers verified
the same runtime configuration and known-origin control receipt before resumption.

The stopped primary instance was restarted. The terminated Spot worker was replaced
by the same instance type in its original availability zone, using the existing
private image and security group. Its retained root volume is mounted at
/mnt/retained-spot2; its home directory and Docker cache are bound to their original
paths. Its retained scratch disk is mounted at /mnt/paired-study-inputs. No retained
volume was reformatted or deleted. A fresh boot volume is separately costed. The
replacement SSH key was checked against authenticated Systems Manager output.

The incomplete deeper comparison directory was moved intact to
/home/ec2-user/paired-study-20261024/interrupted-attempts/halted-SRR36876827_25m.
It will be recomputed from the same verified input rather than mixed with partial
outputs. All already-completed archives remain unchanged. Execution metadata records
the replacement instance and the previous placement.

The new 24-hour guard bounds additional compute at $144 using the two conservative
rate ceilings. Including the prior $47.94 bound gives $191.94, below the existing
$200 working compute allocation. This guard is a cost safeguard, not a user deadline.
The existing account quota is 48 standard on-demand vCPUs and 128 standard Spot vCPUs;
this resumption uses 48 of each. No quota increase was requested.

The foreground supervisor commits each accepted unit after private preservation,
records resource-time costs, verifies worker shutdown, and closes each compute
interval. Missing units remain missing; nothing is counted as a zero result.
