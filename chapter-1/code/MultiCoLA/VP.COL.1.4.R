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


cat("~~~~~~~~~~~Variation partitioning~~~~~~~~~~~\n")
cat('Use the function as follows:\n')
cat('        VP.1.taxa<-VP.COL(M,ENV,Type,typem)\n\n')
cat('        M=output from taxa.pooler,\n')
cat('        ENV = sample by 4 env par (=columns) table,\n')
cat('        Type of cutoff? (all dataset-,"ADS", or sample-,"SAM", based)\n')
cat('        typem = "dominant" or "rare"?\n')
cat('--->type of output:\n') 
cat('	    list with: sum truncated matrix+AdjRsquare\n')
cat('	    VP values\n')


VP.COL<-function(MM,ENV,Type,typem="dominant"){
require(vegan)
  OUTP=readline("\nOutput as text files? (y/n)...\t")
  PLOT=readline("\nPlot the results? (y/n)...\t")
                       
	COP<-function(ODS,z,ENV,Type,typem="dominant"){###################
			#vector to store sum,adjRsq 
			res1<-matrix(NA,1,2)

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
				res1[,1]<-sum(Q3)
			} #end"ADS"
			
			##sample-based cutoff
			if(Type=="SAM"){ 	
				Q1<-ODS
				if(typem=="dominant"){Q1[Q1<z]<-0}	# all species presents less than j times =0
				if(typem=="rare"){Q1[Q1>z]<-0}	# all species presents more than j times =0
				Q3<-CLcol(CLrow(Q1))	#remove rows and columns whose sum=0
				res1[1,1]<-sum(Q3)
			} #end "SAM"                        
			                
##########################
			if (res1[,1]==0){
				res1[,2]<-NA
				VP_Rsq<-matrix(NA,nrow(ODS),1)
				List=list(res1,VP_Rsq,Q3)
				names(List)=c("res1","VP_Rsq","cutoff.table")
			} #end if
      # to avoid conflicts comparing original dataset/NA
			else {   ###else #1
        if (length(Q3)<=nrow(ODS)){
					res1[,2]<-NA
					VP_Rsq<-matrix(NA,nrow(ODS),1)
					List=list(res1,VP_Rsq,Q3)
					names(List)=c("res1","VP_Rsq","cutoff.table")
			 } #end if 
			else {    ###else #2
        if (nrow(Q3)<nrow(ENV)) {
					res1[,2]<-NA
					VP_Rsq<-matrix(NA,nrow(ODS),1)
					List=list(res1,VP_Rsq,Q3)
					names(List)=c("res1","VP_Rsq","cutoff.table")
			 }
				else { ###else #3
					Q2<-decostand(Q3,"hel")[1:nrow(Q3),1:ncol(Q3)]
					#####################################
					#Variation partitioning
					###Transform input as model matrices
					ENV1<-model.matrix(~.,as.data.frame(ENV[,1]))[,-1]
					ENV2<-model.matrix(~.,as.data.frame(ENV[,2]))[,-1]
					ENV3<-model.matrix(~.,as.data.frame(ENV[,3]))[,-1]
					ENV4<-model.matrix(~.,as.data.frame(ENV[,4]))[,-1]
					Q2mod<-model.matrix(~.,as.data.frame(Q2))[,-1]
	
					###Variation partitioning with 4 variables 
					VP_Q2mod<-varpart(Q2mod,ENV1,ENV2,ENV3,ENV4)
					VP_Rsq<-as.data.frame(VP_Q2mod$part$indfract$Adj.R.square)
					res1[,2]<-VP_Q2mod$part$fract[c("[abcdefghijklmno] = All"),c("Adj.R.square")]
					List<-list(res1,VP_Rsq,Q3)
					names(List)<-c("res1","VP_Rsq","cutoff.table")
				} #end else #3
			} #end else #2
     } #end else #1		
return(List)
} #end COP


	if(Type=="ADS"){
	VP.taxa<-function(MM,ENV,Type){
		VPcutoff<-function(ODS,ENV,Type){
		  #create a matrix to store VPvalues for each CO
		  result1<-matrix(NA,16,21)	
		  a<-colnames(ENV)[1]
		  b<-colnames(ENV)[2]
		  d<-colnames(ENV)[3]
		  e<-colnames(ENV)[4]
		  row.names(result1)<-c(a,b,d,e,paste(a,b,sep="+"),paste(a,d,sep="+"),paste(a,e,sep="+"),paste(b,d,sep="+"),paste(b,e,sep="+"),paste(d,e,sep="+"),paste(a,b,d,sep="+"),paste(a,b,e,sep="+"),paste(a,d,e,sep="+"),paste(b,d,e,sep="+"),"All","Unexplained")
		  colnames(result1)<-c(0.01,seq(0.05,0.95,by=0.05),0.99)
		  #create a matrix to store sum,adjRsq,envAIC of all cutoffs
		  result2<-matrix(NA,21,2)	
		  colnames(result2)<-c("Sum","Adj.R.square")
		  row.names(result2)<-c(0.01,seq(0.05,0.95,by=0.05),0.99)
		  LISTRES<-vector("list",21)
		  names(LISTRES)=c(0.01,seq(0.05,0.95,by=0.05),0.99)
		result3<-vector("list",21)
			names(result3)=c(0.01,seq(0.05,0.95,by=0.05),0.99)
			for(i in 1:21){
       LISTRES[[i]]=COP(ODS,z=as.numeric(names(LISTRES)[i]),ENV,Type,typem)   
			 result1[,(22-i)]<-LISTRES[[c(i,2)]][1:16,]
			 result1[,i][result1[,i]<0]<-0
			 result2[(22-i),]<-LISTRES[[c(i,1)]]
			 result3[[i]]<-LISTRES[[c(i,3)]]
		  }
	  LIST2<-list(result1,result2,result3)
	  names(LIST2)<-c("VP_Rsq","res1","cutoff.tables")
		return(LIST2)
		} #end VPcutoff
		
  list.ecol<-vector("list",length(MM))
  names(list.ecol)<-names(MM) 	
  for(j in 1:length(MM)){		
    list.ecol[[j]]<-VPcutoff(MM[[j]],ENV,Type)
  } #end for
  return(list.ecol)
	} #end VP.taxa
} #end if "ADS"
	
	if(Type=="SAM"){
	VP.taxa<-function(MM,ENV,Type){
    limSAMco=as.numeric(readline("\nIf SAM-based only, maximum cutoff value? (e.g. 208)...\t"))
		VPcutoff<-function(ODS,ENV,Type){
      SAM_perc<-limSAMco*c(0.005,0.01,0.015,0.025,0.05,0.075,0.1,0.15,0.25,0.4,0.5,0.6,0.75,0.85,1)
		  #create a matrix to store VPvalues for each CO
		  result1<-matrix(NA,16,length(SAM_perc))	
		  a<-colnames(ENV)[1]
		  b<-colnames(ENV)[2]
		  d<-colnames(ENV)[3]
		  e<-colnames(ENV)[4]
		  row.names(result1)<-c(a,b,d,e,paste(a,b,sep="+"),paste(a,d,sep="+"),paste(a,e,sep="+"),paste(b,d,sep="+"),paste(b,e,sep="+"),paste(d,e,sep="+"),paste(a,b,d,sep="+"),paste(a,b,e,sep="+"),paste(a,d,e,sep="+"),paste(b,d,e,sep="+"),"All","Unexplained")
      colnames(result1)<-round(SAM_perc,0)
		  #create a matrix to store sum,adjRsq,envAIC of all cutoffs
		  result2<-matrix(NA,length(SAM_perc),2)	
		  colnames(result2)<-c("Sum","Adj.R.square")
      row.names(result2)<-round(SAM_perc,0)
		  LISTRES<-vector("list",length(SAM_perc))
		  names(LISTRES)<-SAM_perc
		result3<-vector("list",15)
			names(result3)<-round(SAM_perc,0)
			for(i in 1:length(SAM_perc)){
       LISTRES[[i]]=COP(ODS,z=as.numeric(names(LISTRES)[i]),ENV,Type,typem)   
			 result1[,i]<-LISTRES[[c(i,2)]][1:nrow(ODS),]
			 result1[,i][result1[,i]<0]<-0
			 result2[i,]<-LISTRES[[c(i,1)]]
			result3[[i]]<-LISTRES[[c(i,3)]]
		  } #end for
		  
	  LIST2<-list(result1,result2,result3)
	  names(LIST2)<-c("VP_Rsq","res1","cutoff.tables")
		return(LIST2)
		} #end VPcutoff
		
  list.ecol<-vector("list",length(MM))
  names(list.ecol)<-names(MM) 	
  for(j in 1:length(MM)){		
    list.ecol[[j]]<-VPcutoff(MM[[j]],ENV,Type)
  } #end for                 
  return(list.ecol)
	} #end VP.taxa
} #end if "SAM"
        
result.VP<-VP.taxa(MM,ENV,Type)

  if(OUTP=="y"){
      for(i in 1:length(MM)){
        write.table(result.VP[[c(i,2)]],paste(names(MM[i]),".sum.adjRsq.txt",sep=""),quote=FALSE)
        write.table(result.VP[[c(i,1)]],paste(names(MM[i]),".VarPart.txt",sep=""),quote=FALSE)
      }
  } #end if
  
  if(PLOT=="y"){
    par(mfrow=c(round(length(MM)/2,0),2),oma = c(0, 1, 2, 0),mar=c(5,4,4,2))
    
    for(i in 1:length(MM)){
      barplot(result.VP[[c(i,1)]][1:15,]*100,ylim=c(0,100),col=c("red","darkorange","maroon","royalblue","purple","lightgrey","orchid","green","lightpink","lightblue","brown","yellow","seagreen","blue","tomato")
,xlab=c("%cutoff removed"),ylab=c("% of variation explained"),border=NA)
	title(paste(names(MM[i]),sep=""))
    }
    #with(1, legend(0.1,1,row.names(result.VP[[c(1,1)]][1:15,]),fill=c("red","darkorange","maroon","royalblue","purple","lightgrey","orchid","green","lightpink","lightblue","brown","yellow","seagreen","blue","tomato"),y.intersp=0.7))
    plot(1, type="n", bty="n",axes=F, xlab="", ylab="")
	legend("bottomright",row.names(result.VP[[c(1,1)]][1:15,]),fill=c("red","darkorange","maroon","royalblue","purple","lightgrey","orchid","green","lightpink","lightblue","brown","yellow","seagreen","blue","tomato"),y.intersp=0.7,ncol=3,cex=0.7)

  mtext("Variation partitioning at all taxonomic levels", outer = TRUE, cex = 1.5)
  
  } #end if    
          
return(result.VP)

} #end of VP.COL


