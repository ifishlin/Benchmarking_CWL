# CWL table of the benchmarking project.

The location of the Trimming CWL
```
/tools/trimming
```

The location of the workflows CWL

```
/workflows
```

1). BAT
| Workflow     | Protocols       | Trimming CWL | workflow CWL |
| ------------- |-------------| -----:| :----: |
| BAT | WGBS | trim_galore.cwl | BAT_workflow.cwl |
| BAT | Swift | trim_galore.cwl | BAT_workflow.cwl |
| BAT | T-WGBS | trim_galore.cwl | BAT_workflow.cwl |
| BAT | EMseq | trim_galore.cwl | BAT_workflow.cwl |

2). Biscuit
| Workflow     | Protocols       | Trimming | workflow CWL |
| ------------- |-------------| -----:| :----: |
| Biscuit | WGBS | trim_galore.cwl | biscuit_singlelib.cwl |
| Biscuit | Swift | trim_galore.cwl | biscuit_singlelib.cwl |
| Biscuit | T-WGBS | trim_galore.cwl | biscuit_multilib.cwl |
| Biscuit | EMseq | trim_galore.cwl | biscuit_singlelib.cwl |
| Biscuit | PBAT | trim_galore.cwl | biscuit_singlelib.cwl |

3). Bismark
| Workflow     | Protocols       | Trimming | workflow CWL |
| ------------- |-------------| -----:| :----: |
| Bismark | WGBS | trim_galore.cwl | bismark_singlelib/CWL/workflows/Bismark_start_with_trimmed.cwl |
| Bismark | Swift | trim_galore.cwl | bismark_singlelib/CWL/workflows/Bismark_start_with_trimmed.cwl |
| Bismark | T-WGBS | trim_galore.cwl | bismark_multilib/CWL/workflows/Bismark_TWGBS_multilib.cwl |
| Bismark | EMseq | trim_galore.cwl | bismark_singlelib/CWL/workflows/Bismark_start_with_trimmed.cwl |
| Bismark | PBAT | trim_galore.cwl | bismark_singlelib/CWL/workflows/Bismark_start_with_trimmed.cwl |

4). BSBolt
| Workflow     | Protocols       | Trimming CWL | workflow CWL |
| ------------- |-------------| -----:| :----: |
| BSBolt | WGBS | trim_galore.cwl | bsbolt_singlelib.cwl |
| BSBolt | Swift | trim_galore.cwl | bsbolt_singlelib.cwl |
| BSBolt | T-WGBS | trim_galore.cwl | bsbolt_multilib.cwl |
| BSBolt | EMseq | trim_galore.cwl | bsbolt_singlelib.cwl |
| BSBolt | PBAT | trim_galore.cwl | bsbolt_singlelib.cwl |

5). bwa-meth
| Workflow     | Protocols       | Trimming CWL | workflow CWL |
| ------------- |-------------| -----:| :----: |
| bwa-meth | WGBS | trimmomatic.cwl | bwameth_singlelib.cwl |
| bwa-meth | Swift | trimmomatic_Swift.cwl | bwameth_singlelib.cwl |
| bwa-meth | T-WGBS | trimmomatic.cwl | bwameth_multilib.cwl |
| bwa-meth | EMseq | trimmomatic.cwl | bwameth_singlelib.cwl |
| bwa-meth | PBAT | trimmomatic_PBAT.cwl | bwameth_singlelib.cwl |

6). FAME
| Workflow     | Protocols       | Trimming CWL | workflow CWL |
| ------------- |-------------| -----:| :----: |
| FAME | WGBS | trim_galore.cwl | fame_workflow.cwl |
| FAME | Swift | trim_galore.cwl | fame_workflow.cwl |
| FAME | T-WGBS | trim_galore.cwl | fame_workflow.cwl |
| FAME | EMseq | trim_galore.cwl | fame_workflow.cwl |
| FAME | PBAT | trim_galore.cwl | fame_workflow.cwl |

7). GSNAP
| Workflow     | Protocols       | Trimming CWL | workflow CWL |
| ------------- |-------------| -----:| :----: |
| GSNAP | WGBS | trim_galore.cwl | gsnap_singlelib.cwl |
| GSNAP | Swift | trim_galore.cwl | gsnap_singlelib.cwl |
| GSNAP | T-WGBS | trim_galore.cwl | gsnap_multilib.cwl |
| GSNAP | EMseq | trim_galore.cwl | gsnap_singlelib.cwl |
| GSNAP | PBAT | trim_galore.cwl | gsnap_singlelib.cwl |

8). methylCtools
| Workflow     | Protocols       | Trimming CWL | workflow CWL |
| ------------- |-------------| -----:| :----: |
| methylCtools | WGBS | trimmomatic.cwl | methylCtools_singlelib.cwl |
| methylCtools | Swift | trimmomatic_Swift.cwl | methylCtools_singlelib.cwl |
| methylCtools| T-WGBS | trimmomatic.cwl | methylCtools_multilib.cwl |
| methylCtools | EMseq | trimmomatic.cwl | methylCtools_singlelib.cwl |
| methylCtools | PBAT | trimmomatic_PBAT_methylCtools.cwl | methylCtools_singlelib.cwl |

9). methylpy
| Workflow     | Protocols       | Trimming CWL | workflow CWL |
| ------------- |-------------| -----:| :----: |
| methylpy | WGBS | trim_galore.cwl | methylpy_singlelib.cwl |
| methylpy | Swift | trim_galore.cwl | methylpy_singlelib.cwl |
| methylpy | T-WGBS | trim_galore.cwl | methylpy_multilib.cwl |
| methylpy | EMseq | trim_galore.cwl | methylpy_singlelib.cwl |
| methylpy | PBAT | trim_galore.cwl | methylpy_singlelib.cwl |
