
mod_my_p_val <- function(val){
  if (is.na(val) == FALSE){
    
    if (val >= 0.05){
      out.p = "p > 0.05"
    } else if (val >= 0.01 && val < 0.05){
      out.p = "p < 0.05"
    } else if (val < 0.01){
      out.p = "p < 0.01"
    }
    
  } else {
    out.p = "NA"
  }
  
  
  return(out.p)
}

type_significant <- function(p_val){
  if (p_val > 0.05){
    out.note = "no significant difference"
  } else if (p_val <= 0.05){
    out.note = "a significant difference"
  }
  return(out.note)
}

significant_cor <- function(p_val){
  if (is.na(p_val) == FALSE){
    
  if (p_val > 0.05){
    out.note = "no significant correlation"
  } else if (p_val <= 0.05){
    out.note = "a significant correlation"
  }

  } else {
    out.note = "the model fitting failed"
    
  }
  return(out.note)
}

significant_mod <- function(p_val){
  if (is.na(p_val) == FALSE){
    
    if (p_val > 0.05){
      out.note = "the model did not significantly explain variation in community composition"
    } else if (p_val <= 0.05){
      out.note = "the model significantly explained variation in community composition"
    }
  } else {
    out.note = "the model fitting failed"
  }
  
  return(out.note)
}



lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


#Function for mean and SD per pathway per group (taken from http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


how.rich <- function(phylo){
  
  comm = data.frame(t(otu_table(phylo)), check.names=FALSE)
  tree = phy_tree(phylo)
  alpha.df = pd(comm, tree, include.root = TRUE)
  
  return(alpha.df)
}

#phylo = pro
get_my_long_lat <- function(phylo){
  
  #get metadata
  meta <- data.frame(sample_data(phylo))
  
  #put a negative sign in front of all the south latitude coords
  for (i in 1:nrow(meta)){
    if (meta$Pole[i] == "Antarctic"){
      meta$Latitude[i] <- c(paste0("-", meta$Latitude[i]))
    }
  }
  
  #Some of our cryoconite holes don't have specific coords so we'll use the glacier coordinates
  #remove letters from glacier coord strings
  # meta$GlacierLatitude <- as.numeric(str_sub(meta$GlacierLatitude,1,nchar(meta$GlacierLatitude)-1))
  # meta$GlacierLongitude <- as.numeric(str_sub(meta$GlacierLongitude,1,nchar(meta$GlacierLongitude)-1))
  
  # glac_lat <- meta$Latitude
  # glac_long <- meta$Longitude
  
  #get latitude of each sample
  Latitude <- as.numeric(meta$Latitude)
  #get longitude of each sample
  Longitude <- as.numeric(meta$Longitude)
  
  #replace any missing cryoconite coords with those of it's glacier
  # for (j in 1:length(Latitude)){
  #   if (is.na(Latitude[j])){
  #     Latitude[j] <- glac_lat[j]
  #     Longitude[j] <- glac_long[j]
  #   }
  # }
  
  #put into matrix
  long.lat <- as.data.frame(cbind(Longitude, Latitude))
  rownames(long.lat) = rownames(meta)
  
  return(long.lat)
  
}


#phylo = pro.ant
meta_df <- function(phylo){
  
  #Get metadata of samples
  samp_data <- data.frame(sample_data(phylo))
  
  #list of variables we're interested in 
  keeps1 <- c("WaterDepth", "SedimentDepth", "TotalDepth", "pH", "DOC", "SO4","Ca", "Area")
  
  keeps2 <-c("DistanceToSea", "Elevation")
  
  samp_data1 <- samp_data[ , (names(samp_data) %in% keeps1)]
  samp_data2 <- samp_data[ , (names(samp_data) %in% keeps2)]
  
  #remove columns with more than 50% missing variables
  samp_data_trim1 <- samp_data1[ lapply( samp_data1, function(x) sum(is.na(x)) / length(x) ) < 0.5 ]
  samp_data_trim2 <- samp_data2[ lapply( samp_data2, function(x) sum(is.na(x)) / length(x) ) < 0.7 ]
  
  #Only keep complete cases
  samp_data_trim_complete1 <- samp_data_trim1[complete.cases(samp_data_trim1),]
  samp_data_trim_complete2 <- samp_data_trim2[complete.cases(samp_data_trim2),]
  
  data_to_test = list(samp_data_trim_complete1, samp_data_trim_complete2)
  
  return(data_to_test)
}


meta_df_combi <- function(phylo){
  
  #Get metadata of samples
  samp_data <- data.frame(sample_data(phylo))
  
  #list of variables we're interested in 
  keeps <- c("CryoconiteLatitude", "TotalDepth", "WaterDepth", "SedimentDepth", "pH", "DOC_mg",
             "SO4_mg", "Ca_mg", "Area", "DistanceToSea", "Elevation")
  
  x <- samp_data[ , (names(samp_data) %in% keeps)]
  
  #get cryoconite hole areas
  # area1 <- pi*samp_data$Radius^2
  # area2 <- pi*(samp_data$NS/2)*(samp_data$EW/2)
  # samp_data$Area <- coalesce(area1,area2)
  # 
  # drops <- c("EW", "NS", "Radius")
  # samp_data <- samp_data[ , !(names(samp_data) %in% drops)]
  
  #make sure data are numeric
  #samp_data[1:17] <- data.frame(lapply(samp_data[1:17],as.numeric))
  
  #change environ names 
  names(x) <- gsub("_", "", names(x))
  names(x) <- gsub("mg", "", names(x))
  names(x) <- gsub(" ", "", names(x))
  
  #change environ names 
  #names(samp_data) <- c("Latitude", "Distance To Sea","Elevation", "Water Depth", "Sediment Depth", "Total Depth","Conductivity", "pH", "DOC", "Cl", "SO4", "Na", "K", "Mg","Ca", "HC03", "Area")
  
  
  return(x)
}

agg_my_counts <- function(phylo, level){
  
  spe <- data.frame(t(otu_table(phylo)))
  
  #get out taxonomy table
  tax_tab <- data.frame(tax_table(phylo))
  #replace NAs with Unknown
  tax_tab[is.na(tax_tab)] <- "Unknown"
  colnames(spe) <- tax_tab[, level]
  
  #aggregate count data
  x <- t(spe)
  new_df=aggregate(x, by=list(rownames(x)),sum)
  new_counts <- t(new_df)
  colnames(new_counts) <- new_counts[1,]
  new_counts <- new_counts[-1,]
  df <- data.frame(new_counts)
  #make dataframe numeric
  df[] <- lapply(df, as.numeric)
  df <- df[,! names(df) %in% c("Unknown")]
  
  #remove rows of the dataframe with zero counts
  df = df[rowSums(df[])>0,]
  
  
  return(df)
}


merge_my_phylo <- function(phylo1, phylo2){
  
  #Extract metadata
  meta1 = data.frame(sample_data(phylo1))
  meta2 = data.frame(sample_data(phylo2))
  
  #Extract ASV table of phylo 1 and 2
  tab1 = data.frame(otu_table(phylo1), check.names = FALSE)
  tab2 = data.frame(otu_table(phylo2), check.names=FALSE)
  
  #change sample names to match
  names(tab1) = meta1$Name
  names(tab2) = meta2$Name
  
  #Where samples are missing, replace with zero counts
  in1not2 = setdiff(names(tab1), names(tab2))
  in2not1 = setdiff(names(tab2), names(tab1))
  df1 =  data.frame(matrix(0, ncol = length(in1not2), nrow = nrow(tab2)))
  df2 =  data.frame(matrix(0, ncol = length(in2not1), nrow = nrow(tab1)))
  names(df1) = in1not2
  names(df2) = in2not1
  tab1 = cbind(tab1, df2)
  tab2 = cbind(tab2, df1)
  setdiff(names(tab1), names(tab2))
  setdiff(names(tab2), names(tab1))
  
  #merge ASV tables
  new.df = rbind(tab1, tab2)
  
  #combine tax tables
  new.tax = rbind(data.frame(tax_table(phylo1)), data.frame(tax_table(phylo2)))
  new.tax.mat = as.matrix(new.tax)
  
  #change metadata row names to match samples
  meta.comb = rbind(meta1, meta2)
  row.names(meta.comb) = meta.comb$Name
  
  #make new phyloseq object
  ASV = otu_table(new.df, taxa_are_rows = TRUE)
  TAX = tax_table(new.tax.mat)
  META = sample_data(meta.comb)
  new.phylo = phyloseq(ASV, TAX, META)
  
  return(new.phylo)
  
}

merge_my_phylo_same_pole <- function(phylo1, phylo2){
  
  #Extract metadata
  meta1 = data.frame(sample_data(phylo1))
  meta2 = data.frame(sample_data(phylo2))
  
  #Extract ASV table of phylo 1 and 2
  tab1 = data.frame(otu_table(phylo1), check.names = FALSE)
  tab2 = data.frame(otu_table(phylo2), check.names=FALSE)
  
  #change sample names to match
  names(tab1) = meta1$Name
  names(tab2) = meta2$Name
  
  #Where samples are missing, replace with zero counts
  in1not2 = setdiff(names(tab1), names(tab2))
  in2not1 = setdiff(names(tab2), names(tab1))
  df1 =  data.frame(matrix(0, ncol = length(in1not2), nrow = nrow(tab2)))
  df2 =  data.frame(matrix(0, ncol = length(in2not1), nrow = nrow(tab1)))
  names(df1) = in1not2
  names(df2) = in2not1
  tab1 = cbind(tab1, df2)
  tab2 = cbind(tab2, df1)
  setdiff(names(tab1), names(tab2))
  setdiff(names(tab2), names(tab1))
  
  #merge ASV tables
  new.df = rbind(tab1, tab2)
  
  #combine tax tables
  new.tax = rbind(data.frame(tax_table(phylo1)), data.frame(tax_table(phylo2)))
  new.tax.mat = as.matrix(new.tax)
  
  #change metadata row names to match samples
  meta2a = meta2[meta2$Name %in% setdiff(meta2$Name, meta1$Name), ]
  meta.comb = dplyr::bind_rows(meta1, meta2a)
  row.names(meta.comb) = meta.comb$Name
  
  # length(unique(c(meta1$SampleID, meta2$SampleID)))
  # length(unique(c(meta1$Name, meta2$Name)))
  
  #make new phyloseq object
  ASV = otu_table(new.df, taxa_are_rows = TRUE)
  TAX = tax_table(new.tax.mat)
  META = sample_data(meta.comb)
  new.phylo = phyloseq(ASV, TAX, META)
  
  return(new.phylo)
  
}
