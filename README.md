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
remotes::install_github("SorenHeidelbach/ModFind")
```


## Example usage

```r
# required inputs
nat_mapping = "nat/mappings.sort.bam"
pcr_mapping = "pcr/mappings.sort.bam"
nat_signal = "nat/signal_mappings.hdf5"
pcr_signal = "pcr/signal_mappings.hdf5"
chunk_size = 1e5

# HDF5 list
h5_list <- list(
  nat = hdf5r::H5File$new(nat_signal, mode = "r"),
  pcr = hdf5r::H5File$new(pcr_signal, mode = "r")
)

# Preprocess
metainfo <- prepare_metainfo(
  nat_mapping,
  pcr_mapping,
  hdf5 = h5_list,
  chunk_size = chunk_size
)

# process the first chunk
chunk <- process_chunk(
    metainfo[[1]],
    h5_list = h5_list,
    chunk_size = chunk_size
  )
```