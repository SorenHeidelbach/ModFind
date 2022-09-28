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

# Process chunked reads
preprocess_all_chunks(
  metainfo,
  h5_list,
  out
)

# Load processed chunks
chunks <- load_processed_chunks(out)

# Embed and cluster processed chunks
clusters <- embed_and_cluster(
  chunks,
  ref_path = reference_path,
  out = out,
  umap_args = list(n_components = 3),
  hdbscan_args = list(minPts = 10)
)

# Get cluster events sequences
setorder(clusters, HDBSCAN)
cluster_sequences <- lapply(
    split(clusters[!is.na(HDBSCAN)], by = "HDBSCAN"),
    function(x) {
      toupper(x$seq)
    }
  )

# Calculate entropy and select motifs
lapply(
    cluster_sequences,
    function(x) calculate_bit_score(x, min_entropy = 1)
  )

# Visualise embedd + clustering
viz <- visualise_clusters(clusters)
ggsave(
  paste_path(out, "cluster_scatter.png"),
  width = 6,
  height = 5,
  viz[[1]]
)
ggsave(
  paste_path(out, "cluster_logo.png"),
  width = 6,
  height = 5,
  viz[[2]]
)

```
