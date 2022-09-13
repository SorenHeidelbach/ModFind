# ModFind

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

## Usage

There are three main functionalities:
- Preprocess; which reformats and gathers mappings of NAT and PCR reads
- Current difference; which calculate statistics between NAT and PCR reads at each reference position and selects candidate modified positions
- Identify motifs; which embeds and cluster candidate positions to identify sequence motifs




