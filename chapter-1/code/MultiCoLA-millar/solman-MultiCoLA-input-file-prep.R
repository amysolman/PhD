MultiCoLA_input_file_prep <- function(my_phylo){
  #extract ASV count table
  counts <- data.frame(otu_table(my_phylo), check.names = FALSE)
  #extract taxonomy
  taxa <- data.frame(tax_table(my_phylo))
  #remove kingdom column
  taxa <- subset(taxa, select=-c(Kingdom))
  #merge and return dataframe
  counts_taxa <- cbind(counts, taxa)
  
  return(counts_taxa)
}