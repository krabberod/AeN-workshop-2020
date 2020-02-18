#' @title Estimate phylogenetic diversity
#'
#' @description This function estimates the Faiths phylogenetic diverstiy from an OTU table and phylogenetic tree
#' @param phylo_seq_object A phyloseq object with an OTU table and phylogenetic tree.
#' @details
#' Faiths PD-index is supported.
#' @return A data.frame of phylogenetic diversity metrics.
#' @keywords phyloseq
#' @export
#' @examples
#' estimate_pd()

estimate_pd <- function(phylo_seq_object){
  OTU <- phyloseq::otu_table(phylo_seq_object)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }

  otutable <- as(OTU, "matrix")
  tree <- phyloseq::phy_tree(phylo_seq_object)

  # Print status message
  message("Calculating Faiths PD-index...")

  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }

  pdtable <- picante::pd(otutable, tree, include.root = FALSE)

  return(pdtable)
}
