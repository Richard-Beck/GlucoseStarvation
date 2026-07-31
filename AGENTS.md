Read docs/project_map.txt when you need to understand the project and the organisation of key files.
project_map.txt is intended as a living document, so if you make any changes to the codebase, determine whether this necessitates an update to project_map.txt. project_map.txt is intended to be relatively terse, so do not bloat it with uneccessary information.
If your current environment is /share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/GlucoseStarvation you are in the HPC system! HPC system suggestions:
1) Run all R scripts via scripts/agentRrunner.sh myScript arg1 arg2... 
2) Run small jobs/scripts directly from the terminal.
3) For larger jobs/scripts (parallel/batch, requiring GPU, large RAM, expected execution time > 5 mins) use SLURM.

When using SLURM, set QOS strategically to ensure minimum expected execution time for the job. See:
Current QOS settings from
`sacctmgr show qos format=Name,Priority,MaxTRESPU,MaxJobsPU,MaxSubmitJobsPU,GrpTRES -P`:

```text
Name|Priority|MaxTRESPU|MaxJobsPU|MaxSubmitPU|GrpTRES
normal|0|cpu=64|1000001||
small|0|cpu=475|1000001||
large|0|cpu=1425|1000001||
xlarge|0|cpu=1800|1000001||
partsmall|0||1000001||cpu=250
medium|0|cpu=950|1000001||
xxlarge|0|cpu=3000|1000001||
```

