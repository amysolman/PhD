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

COtables<-function(ODS,Type="ADS",typem="dominant"){                                       
	COP<-function(ODS,z,Type,typem){###################
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
				Q3<-Q1[,-which(apply(Q1,2,function(x)all(is.na(x))))]
				Q3<-CLcol(CLrow(as.data.frame(Q3)))	#remove rows and columns whose sum=0
			} #end"ADS"
			
			##sample-based cutoff
			if(Type=="SAM"){ 	
				Q1<-ODS
				if(typem=="dominant"){Q1[Q1<z]<-0}	# all species presents less than j times =0
				if(typem=="rare"){Q1[Q1>z]<-0}	# all species presents more than j times =0
				Q3<-CLcol(CLrow(Q1))	#remove rows and columns whose sum=0
			} #end "SAM"                        
return(Q3)
} #end COP

	if(Type=="ADS"){
		  #create a matrix to store VPvalues for each CO
		  LISTRES<-vector("list",21)
      names(LISTRES)<-c(0.01,seq(0.05,0.95,by=0.05),0.99)
			for(i in 1:21){
       LISTRES[[i]]=COP(ODS,z=as.numeric(names(LISTRES)[i]),Type,typem)   
		  }
	} #end if "ADS"
	
	if(Type=="SAM"){
      #create a matrix to store VPvalues for each CO
		  LISTRES<-vector("list",15)
      names(LISTRES)<-c(1,2,3,5,10,15,20,30,55,80,105,130,155,180,208)
			for(i in 1:15){
       LISTRES[[i]]=COP(ODS,z=as.numeric(names(LISTRES)[i]),Type,typem)   
		  }
	} #end if "SAM"
	return(LISTRES)
} #end of VP.COL
###################################################################################



			
