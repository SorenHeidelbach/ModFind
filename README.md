# nanomotif

## Overview

Theese functions process the output of megalodon and identify modified motifs. 
It depends on signal and read mappings from both natural reads (NAT) and non-modified reads (PCR).

File inputs are 

- signal_mappings to reference/assembly (NAT and PCR)
- read_mappings to reference/assembly (NAT and PCR)
- Reference/assembly sequence

## Installing

```r
install.packages("remotes")
remotes::install_github("SorenHeidelbach/nanomotif")
```


## Example usage

```r
# required inputs
nat_mapping = "nat/mappings.sort.bam"
pcr_mapping = "pcr/mappings.sort.bam"
nat_signal = "nat/signal_mappings.hdf5"
pcr_signal = "pcr/signal_mappings.hdf5"
reference_path = "path/to/reference.fasta"
out = "path/to/output"

nanomotif(
    nat_mapping = nat_mapping,
    pcr_mapping = pcr_mapping,
    nat_signal = nat_signal,
    pcr_signal = pcr_signal,
    reference_path = reference_path
    out = out,
)


```
