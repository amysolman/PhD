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

cat("~~~~~~~~~~~Abundance, Correlation and Procrustes~~~~~~~~~~~\n")
cat("~~~~~~~~~~~To obtain a suitable output for figures~~~~~~~~~~~\n")
cat('Use the function as follows:\n')
cat('        output.all<-cutoff.impact.fig(M)\n')
cat('        M=output from ecology.extractor,\n')
cat('--->type of output:\n')
cat('	    vector with: total sum of pyro-tags, correlation coefficient,\n')
cat('	    procrustes R value\n')

cutoff.impact.fig<-function(M){
  OUTP=readline("\nOutput as text files? (y/n)...\t")
  PLOT=readline("\nPlot the results? (y/n)...\t")
  list.all<-vector("list",3)
  names(list.all)<-c("Abundance","Non-par.correlation","Procrustes")
  list.all[[1]]<-matrix(NA,nrow(M[[1]]),length(M))
  list.all[[2]]<-matrix(NA,nrow(M[[1]]),length(M))
  list.all[[3]]<-matrix(NA,nrow(M[[1]]),length(M))
  colnames(list.all[[1]])<-names(M)
  colnames(list.all[[2]])<-names(M)
  colnames(list.all[[3]])<-names(M)
  row.names(list.all[[1]])<-row.names(M[[1]])
  row.names(list.all[[2]])<-row.names(M[[1]])
  row.names(list.all[[3]])<-row.names(M[[1]])
  
  for(i in 1:length(M)){
    list.all[[1]][,i]<-cbind(M[[i]][,1])
    list.all[[2]][,i]<-cbind(M[[i]][,2])
    list.all[[3]][,i]<-cbind(M[[i]][,3])
  } #end for

  if(OUTP=="y"){
  write.table(list.all[[1]],"abundance.txt",quote=FALSE)
  write.table(list.all[[2]],"non-par.correlation.txt",quote=FALSE)
  write.table(list.all[[3]],"procrustes.txt",quote=FALSE)
  } #end if
  
  if(PLOT=="y"){
    par(mfrow=c(3,1))
    plot(row.names(list.all[[1]]),list.all[[1]][,1],type="l",xlab=c("%cutoff removed"),ylab=c("Abundance in each matrix"))
    for(i in 1:length(M)){
      lines(row.names(list.all[[1]]),list.all[[1]][,i],col=i)
    }
    plot(row.names(list.all[[2]]),list.all[[2]][,1],type="l",ylim=c(0,1),xlab=c("%cutoff removed"),ylab=c("Non-par.correlation"))
    for(i in 1:length(M)){
      lines(row.names(list.all[[2]]),list.all[[2]][,i],col=i)
    }
    plot(row.names(list.all[[3]]),list.all[[3]][,1],type="l",ylim=c(0,1),xlab=c("%cutoff removed"),ylab=c("Procrustes correlation"))
    for(i in 1:length(M)){
      lines(row.names(list.all[[3]]),list.all[[3]][,i],col=i)
      legend(0,0.8,colnames(as.data.frame(list.all[[1]])),col=seq(1:length(M)),lty=1,y.intersp=0.7)
    }
  } #end if
  
return(list.all)
} #end cutoff.impact.fig
