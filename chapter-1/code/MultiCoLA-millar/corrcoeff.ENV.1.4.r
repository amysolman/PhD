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


corrcoeff<-function(SPE,ENV){
require(vegan)
res2<-matrix(NA,1,5)
colnames(res2)<-c("Sum",paste("RDA1.",colnames(ENV),sep=""))
	res2[,1]<-sum(SPE)
	if (nrow(SPE)<nrow(SPE)){
			res2[,2:5]<-c(NA,NA,NA,NA)
		} # to avoid conflicts comparing original dataset/NA
		else {
		ENV1<-model.matrix(~.,as.data.frame(ENV))[,-1]
		Q2<-decostand(SPE,"hel")[1:nrow(SPE),1:ncol(SPE)]
		#RDA1 axis values for each env par
		R1=rda(Q2~ENV1)
		res2[,2:5]=summary(R1)$biplot[,"RDA1"]		
		}	
return(res2)
}
