# CWL implemetation

We implemented 10 methylation calling workflow in Common Workflow Language.

1. BAT
2. Biscuit
3. Bismark
4. BSBolt
5. bwa-meth
6. FAME
7. gemBS
8. GSNAP
9. methylCtools
10. methylpy

## Usage
Each workflow has two version, one for single library and the other for mulitple library readsets. 

You can find two CWLs under the directory with names 
```
<workflow>.singlelib.cwl 

<workflow>.mulitlib.cwl
```

Cwltool and Docker are requried.

ALL workflows utilize docker and you can execute them without any software installation.
Usage:
```
cwltool <workflow>.singlelib.cwl <workflow>.singlelib.yml
```

# The Benchmarking Project

Analysis of bisulfite sequencing data relies upon processing that generally includes four core steps: 
read preprocessing, alignment, post-alignment processing and calling of methylation states. 
An impressive number of tools for each of the steps or their combinations, workflows integrating them as well as turn-key solutions have been proposed.

Despite of this versatility, so far only few attempts have been made to systematically evaluate complete processing workflows in a standardized and unbiased analysis. Previous benchmarks either focused upon a single processing task, e.g. predominantly alignment software (Kunde-Ramamoorthy, Coarfa et al. 2014, Tran, Porter et al. 2014, Sun, Han et al. 2018, Grehl, Wagner et al. 2020). None of the benchmarks covered a substantial number of tools. Most importantly, none of the studies was based on a reasonable gold-standard data set.

To bridge this gap, we set out to perform a thorough benchmarking study of bisulfite sequencing workflows. At the core of our benchmark is a set of samples with highly accurate methylation calls (Bock, Halbritter et al. 2016), which we use as the gold-standard. We evaluate the software in the context of five most widely used sequencing protocols and propose protocol-specific choice of workflows. To simplify the choice of workflows and enable continuity we developed rich data presentation and benchmarking resources. To our knowledge, this is the most comprehensive benchmarking study of bisulfite sequencing to date.
