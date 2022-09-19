

dt_read_mapping <- function(chunk_size = 1e3, type = "nat"){
  dt <- data.table(
    x = as.character(1:1e4),
    type = type,
    strand = c("+", "-")[sample(1:2, 1e4, replace = TRUE)],
    pos = 1:1e4,
    qwidth = sample(50:1e3, 1e3, replace = TRUE)
  )[
      , chunk := pos %/% chunk_size
  ]
}

dt_read_mapping_unequal <- function(chunk_size = 1e3){
  dt_nat <- rbind(
      dt_read_mapping(chunk_size, type = "nat"),
      dt_read_mapping(chunk_size, type = "nat")
    )
  dt_pcr <- dt_read_mapping(chunk_size, type = "pcr")
  dt <- rbind(dt_nat, dt_pcr)
}