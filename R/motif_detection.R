


embed_umap <- function(features,
                       dims = 2,
                       epochs = 500,
                       neighbors = 15,
                       min_dist = 0.1) {
    # Assert all numeric columns
    assert::assert(
        all(sapply(
            features,
            FUN = function(x) is.numeric(x)
        )),
        msg = "Feature dt contain non numeric columns"
    )
    # Assert no NAs
    assert::assert(
        all(sapply(
            features,
            FUN = function(x) all(!is.na((x)))
        )),
        msg = "Feature dt contain NA in some columns, replace with value"
    )

    embedding <- umap::umap(
        features,
        n_components = dims,
        n_epochs = epochs,
        n_neighbors = neighbors,
        min_dist = min_dist
        )
    dt <- data.table(embedding$layout)
    setnames(dt, paste0("V", 1:dims), paste0("UMAP", 1:dims))
}



cluster_emb <- function(embedding,
                        min_pts = 10){
    # Assert all numeric columns
    assert::assert(
        all(sapply(
            embedding,
            FUN = function(x) is.numeric(x)
        )),
        msg = "Embedding contain non numeric columns"
    )
    # Assert no NAs
    assert::assert(
        all(sapply(
            embedding,
            FUN = function(x) all(!is.na((x)))
        )),
        msg = "Embedding contain NA in some columns, replace with value"
    )

    clusters_hdbscan <- dbscan::hdbscan(embedding, minPts = min_pts)
    embedding[
            , HDBSCAN := as.factor(clusters_hdbscan$cluster)
        ][
            , cluster_prob := clusters_hdbscan$membership_prob
        ]
}


add_event_frame_sequence <- function(sequence_list) {
    
}