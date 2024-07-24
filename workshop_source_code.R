

# Source code for the DAPPG workshop

###########################
# WormCat is from Amy Walker's group
# WormCat: an online tool for annotation and visualization of Caenorhabditis elegans genome-scale data
# Amy D Holdorf, Daniel P Higgins, Anne C. Hart, Peter R Boag, Gregory Pazour, Albertha J. M. Walhout, Amy Karol Walker
# GENETICS February 1, 2020 vol. 214 no. 2 279-294;
##########################

.merger_cats <- function(rgs_annotated_cat, annotated_cat, total_annotations_count, total_rgs_count, n) {
  
  merged_cats <- base::merge(rgs_annotated_cat, annotated_cat, by = "Var1", all.x = TRUE)
  colnames(merged_cats) <- c("Category", "RGS", "AC")
  
  # Step 5: Build contingency table for each category in RGS vs AC
  df <- data.frame(Category = character(),
                   RGS = double(),
                   AC = double(),
                   Fold = double(),
                   PValue = double(),
                   stringsAsFactors = FALSE)
  
  fact_character <- levels(merged_cats$Category)[as.numeric(merged_cats$Category)]
  
  for (i in 1:nrow(merged_cats)) {
    if (is.na(merged_cats$RGS[i]) | is.na(merged_cats$AC[i])) {
      pvalue <- NA
    } else {
      stat <- fisher.test(matrix(c(merged_cats$RGS[i], total_rgs_count,
                                   merged_cats$AC[i],  total_annotations_count),
                                 nrow = 2, ncol = 2),
                          alternative = "greater")
      pvalue <- stat$p.value
    }
    
    df[nrow(df) + 1, ] <- list(Category = fact_character[i],
                               RGS = merged_cats$RGS[i],
                               AC = merged_cats$AC[i],
                               Fold = (merged_cats$RGS[i]/sum(merged_cats$RGS))/(merged_cats$AC[i]/sum(merged_cats$AC)),
                               pvalue)
  }
  
  sorted_df <- df[with(df, order(PValue)), ]
  return(sorted_df)
}

.worm_cat_acceptable_pvalues <- function(rgs_fisher_cat) {
  Bonferroni <- NULL
  
  rgs_fisher_cat <- na.omit(rgs_fisher_cat)
  
  rgs_fisher_cat[order(rgs_fisher_cat$PValue), ]
  
  bonferroni <- p.adjust(rgs_fisher_cat$PValue, method = "bonferroni")
  rgs_fisher_cat <- data.frame(rgs_fisher_cat, Bonferroni = bonferroni)
  
  ### Acceptable is to be 0.01 on Bonferroni
  rgs_fisher_cat <- subset(rgs_fisher_cat, Bonferroni < 0.05)
}

run_worm_cat <- function(gene_set, gene_data, tier) {
  worm_cat_tier <- c("WormCategory.1", "WormCategory.2", "WormCategory.3")[tier]
  
  # Step 1: Get the gene set
  rgs <- data.frame(table(gene_data[gene_set, worm_cat_tier]))
  
  # Step 2: Get the annotations
  annotated_cat <- data.frame(table(gene_data[,worm_cat_tier]))
  
  # Step 3: Get the total number of annotations
  total_annotations_count <- sum(annotated_cat$Freq)
  
  # Step 4: Get the total number of RGS
  total_rgs_count <- sum(rgs$Freq)
  
  # Step 5: Merge the two dataframes
  merged_cats <- .merger_cats(rgs, annotated_cat, total_annotations_count, total_rgs_count, tier)
  
  # Step 6: Get the acceptable p-values
  acceptable_pvalues <- .worm_cat_acceptable_pvalues(merged_cats)
  
  return(acceptable_pvalues)
}






