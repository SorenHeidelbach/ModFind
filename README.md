# nanomotif

## Overview

Theese functions process the output of megalodon and identify modified motifs. 
It depends on singal and read mappings from both natural reads (NAT) and non-modified reads (PCR).
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
chunk_size = 1e5

# Load read metainfo
metainfo <- prepare_metainfo(
  nat_mapping = nat_mapping,
  pcr_mapping = pcr_mapping,
  nat_hdf5 = nat_signal,
  pcr_hdf5 = pcr_signal,
  chunk_size = chunk_size
)

# Process chunked reads
preprocess_all_chunks(
  metainfo = metainfo,
  nat_hdf5 = nat_signal,
  pcr_hdf5 = pcr_signal,
  out = out,
  chunk_size = chunk_size
)

# Identify motifs
find_motifs(
  path_chunk_stats = out,
  path_ref = reference,
  out = out
)
plot_motifs(
  cluster_path = file.path(out, "/1/clusters.tsv"),
  plot_out  = file.path(out, "/1"),
  path_ref = reference,
  chunks_path = out
)
```
