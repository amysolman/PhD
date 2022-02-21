# Copyright (C) 2011 Angelique Gobet & Alban Ramette
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

cat("~~~~~~~~~~~Cutoff impact~~~~~~~~~~~\n")
cat('Use the function as follows:\n')
cat('        storing_name<-cutoff.impact(all_taxa_pooled,Type="ADS",corcoef="spearman",typem="dominant")\n\n')
cat('        M=output from taxa.pooler,\n')
cat('        Type of cutoff? (all dataset-,"ADS", or sample-,"SAM", based)\n')
cat('        Correlation? ("spearman", "kendall", "pearson")\n')
cat('        Which matrix type? ("dominant" or "rare"?)\n')
cat('--->type of output:\n') 
cat('	    list with: total sum of pyro-tags, correlation coefficient,\n')
cat('	    procrustes R value for each cutoff and all taxonomic levels\n')

cutoff.impact<-function(MM,Type="ADS",corcoef="spearman",typem="dominant"){
require(MASS)
require(vegan)
details=readline("\nDetails of the NMDS calculations? (y/n)...\t")
#calculation of the cut-off matrices, correlation coefficient, procrustes R value
CoCalc<-function(ODS,z,Type,corcoef){
	res1<-matrix(NA,1,3)	#create a matrix to store mantel and procrustes data
	colnames(res1)<-c("Sum","corrcoeff","R Procrustes")
	   
        #to remove all the lines in the matrix for which the sum of the line is 0
        CLrow<-function(m) {
          #to create a column of 0 in the last column
          m=cbind(m,matrix(0,nrow(m),1))
          m[,ncol(m)]=apply(m,1,sum)
          #to keep only the lines without 0
          mclean=subset(m,m[,ncol(m)]!=0)
          mclean=mclean[,-ncol(mclean)]
          return (mclean)
         }
         ########################
        CLcol<-function(m) {
          #transpose the matrix and apply the same function as before
          #then transpose back
          m=t(m)
          #to create a column of 0 in the last column
          m=cbind(m,matrix(0,nrow(m),1))
          m[,ncol(m)]=apply(m,1,sum)
          #to keep only the lines without 0
          mclean=subset(m,m[,ncol(m)]!=0)
          mclean=mclean[,-ncol(mclean)]
          mclean=t(mclean)
          return (mclean)
         }

			#Application of a percentage cut-off to the original dataset to obtain abundant dataset
			##all dataset-based cutoff
			if(Type=="ADS"){                   
				M<-rbind(ODS,apply(ODS,2,sum))	#add the column sum as a last row of the matrix M
				if(typem=="dominant"){N<-M[,order(M[nrow(M),],decreasing=TRUE)]}	#order the columns of M by their decreasing sum
    		if(typem=="rare"){N<-M[,order(M[nrow(M),])]}  #order the columns of M by their increasing sum
			   	Q<-N[1:(nrow(N)-1),]	#remove the last row (with the sum of the columns)
	   			L<-ncol(Q)
   				K<-nrow(Q)
  				M1<-t(matrix(NA,L))	#create a vector to store sum of successive matrices
  				Q1<-matrix(NA,K,L)	#create a matrix to store new data
  				perc<-z*sum(ODS)
				  for (i in 1:L){	###for #1
				  	M1[,i]<-sum(Q[,1:i]) 
			   		if (M1[,i]<=perc) {Q1[,1:i]=Q[,1:i]} 
				    	row.names(Q1)=row.names(Q)
			     		colnames(Q1)=colnames(Q)
		  			if (M1[,i]>perc) {Q1[,1:i]==0}
		  		}#end for #1
				Q2<-Q1[,-which(apply(Q1,2,function(x)all(is.na(x))))]
				Q2<-CLcol(CLrow(as.data.frame(Q2)))	#remove rows and columns whose sum=0
				res1[1,1]<-sum(Q2)
      } #end "ADS"
			
			##sample-based cutoff
			if(Type=="SAM"){ 	
				Q1<-ODS
				if(typem=="dominant"){Q1[Q1<z]<-0}	# all species presents less than j times =0
				if(typem=="rare"){Q1[Q1>z]<-0}	# all species presents more than j times =0
				Q2<-CLcol(CLrow(Q1))	#remove rows and columns whose sum=0
				res1[1,1]<-sum(Q2)
		  } #end "SAM"          

	if (res1[1,1]==0){res1[1,2:3]<-cbind(NA,NA)} # to avoid conflicts when comparing original dataset to NA
	else {	if (length(Q2)<=nrow(ODS)){res1[1,2:3]<-cbind(NA,NA)} ###else #1
		else { if (nrow(Q2)<nrow(ODS)) {res1[1,2:3]<-cbind(NA,NA)}
			else { ###else #2
#####################################
	#Correlation and Procrustes calculations
		ODSdist<-vegdist(ODS,method="bray")	#distance matrix of the original dataset
		Q2dist<-vegdist(Q2,distance="bray")	#distance matrix of the truncated dataset
		ODSQ2cor<-cor.test(ODSdist,Q2dist,method=corcoef) #correlation between matrices
		ODSdist2<-ODSdist
		ODSdist2[ODSdist2==0]<-10e-20	#replace 0 by 10e-20 for original dataset
		Q2dist2<-Q2dist
		Q2dist2[Q2dist2==0]<-10e-20	#replace 0 by 10e-20 for truncated dataset
		if(details=="y"){
		  ODSNMDS<-isoMDS(ODSdist2)	#NMDS for original dataset
		  Q2NMDS<-isoMDS(Q2dist2)		#NMDS for truncated dataset
		  }
		  else{if(details=="n"){     ###else #3
		  ODSNMDS<-isoMDS(ODSdist2,trace=0)	#NMDS for original dataset
		  Q2NMDS<-isoMDS(Q2dist2,trace=0)		#NMDS for truncated dataset
		  }}  		  
		ODSQ2procrustes<-protest(ODSNMDS,Q2NMDS) #procrustes
		res1[1,2]<-cbind(ODSQ2cor$estimate)
		res1[1,3]<-cbind(ODSQ2procrustes$t0)
		} #end else #2
	} #end else #1
 }
#####################################
return(res1)
}#end function CoCalc

################################################################
#application of the function COP at different all dataset-based cutoffs
	if(Type=="ADS"){   
    ecol.ext.all<-function(MM,Type,corcoef){   
		  allcorr<-function(ODS,Type,corcoef){ 		
        ADS_perc<-c(0.01,seq(0.05,0.95,by=0.05),0.99)
        table_taxa<-matrix(NA,length(ADS_perc),3)
        row.names(table_taxa)<-ADS_perc
        colnames(table_taxa)<-c("Sum","corrcoeff","R Procrustes")
        for(i in 1:length(ADS_perc)){
				  table_taxa[((length(ADS_perc)+1)-i),]<-CoCalc(ODS,z=as.numeric(ADS_perc[i]),Type,corcoef)
		    }
		    #ROW<-row.names(table_taxa)
        #table_taxa<-table_taxa[order(as.numeric(row.names(table_taxa)),decreasing=TRUE),]
        #row.names(table_taxa)<-ROW
		    return(table_taxa)
		  } #end allcorr

  list.ecol<-vector("list",length(MM))
  names(list.ecol)<-names(MM)
            
  for(j in 1:length(MM)){		
    list.ecol[[j]]<-allcorr(MM[[j]],Type,corcoef)
  }
  return(list.ecol)
  } #end ecol.ext.all
 }#end if "ADS"

#application of the function COP at different sample-based cutoffs
	if(Type=="SAM"){
    ecol.ext.all<-function(MM,Type,corcoef){   
      limSAMco=as.numeric(readline("\nIf SAM-based only, maximum cutoff value? (e.g. 208)...\t"))
  		allcorr<-function(ODS,Type,corcoef){  
        SAM_perc<-limSAMco*c(0.005,0.01,0.015,0.025,0.05,0.075,0.1,0.15,0.25,0.4,0.5,0.6,0.75,0.85,1)
        table_taxa<-matrix(NA,length(SAM_perc),3)
        row.names(table_taxa)<-round(SAM_perc,0)
        colnames(table_taxa)<-c("Sum","corrcoeff","R Procrustes")
        for(i in 1:length(SAM_perc)){
				  table_taxa[i,]<-CoCalc(ODS,z=as.numeric(SAM_perc[i]),Type,corcoef)
        }
		    return(table_taxa)
		  } #end allcorr

  list.ecol<-vector("list",length(MM))
  names(list.ecol)<-names(MM)
	
  for(j in 1:length(MM)){		
    list.ecol[[j]]<-allcorr(MM[[j]],Type,corcoef)
  }
  return(list.ecol)
  } #end ecol.ext.all
 }#end if "SAM"

result.allcorr<-ecol.ext.all(MM,Type,corcoef)
return(result.allcorr)
            
} #end cutoff.impact

