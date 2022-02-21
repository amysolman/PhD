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

signif<-function(SPE,ENV){
require(vegan)
result3<-matrix(NA,1,5)
colnames(result3)<-c("Pw","Pe1","Pe2","Pe3","Pe4")
#colnames(result3)<-c("whole.sig","ENV1.sig","ENV2.sig","ENV3.sig","ENV4.sig")

	###Transform input as model matrices
	if (length(SPE)<=nrow(SPE)){
		result3<-c("NA","NA","NA","NA")            
	} # to avoid conflicts comparing original dataset/NA
	else {
		ENV1<-model.matrix(~.,as.data.frame(ENV[,1]))[,-1]
		ENV2<-model.matrix(~.,as.data.frame(ENV[,2]))[,-1]
		ENV3<-model.matrix(~.,as.data.frame(ENV[,3]))[,-1]
		ENV4<-model.matrix(~.,as.data.frame(ENV[,4]))[,-1]
		Q2<-decostand(SPE,"hel")[1:nrow(SPE),1:ncol(SPE)]
		#significance
		###significance of whole model
		whole<-permutest(rda(Q2~ENV1+ENV2+ENV3+ENV4),permutations=1000)
		Sw<-sort(whole$F.perm)
		if(length(Sw[Sw>whole$F.0])==0){
			result3[,1]<-c("<0.001")
		}
		else{
			result3[,1]<-length(Sw[Sw>whole$F.0])/length(whole$F.perm)
		}		
		###significance of each environmental parameter
		ENV1.sig<-permutest(rda(Q2~ENV1+Condition(ENV2)+Condition(ENV3)+Condition(ENV4)),permutation=1000,model="full")
		Se1<-sort(ENV1.sig$F.perm)
		if(length(Se1[Se1>ENV1.sig$F.0])==0){
			result3[,2]<-c("<0.001")
		}
		else{		
			result3[,2]<-length(Se1[Se1>ENV1.sig$F.0])/length(ENV1.sig$F.perm)		
		}
		ENV2.sig<-permutest(rda(Q2~ENV2+Condition(ENV1)+Condition(ENV3)+Condition(ENV4)),permutation=1000,model="full")
		Se2<-sort(ENV2.sig$F.perm)
		if(length(Se2[Se2>ENV2.sig$F.0])==0){
			result3[,3]<-c("<0.001")
		}
		else{	
			result3[,3]<-length(Se2[Se2>ENV2.sig$F.0])/length(ENV2.sig$F.perm)		
		}
		ENV3.sig<-permutest(rda(Q2~ENV3+Condition(ENV2)+Condition(ENV1)+Condition(ENV4)),permutation=1000,model="full")
		Se3<-sort(ENV3.sig$F.perm)
		if(length(Se3[Se3>ENV3.sig$F.0])==0){
			result3[,4]<-c("<0.001")
		}
		else{	
			result3[,4]<-length(Se3[Se3>ENV3.sig$F.0])/length(ENV3.sig$F.perm)		
		}	
		ENV4.sig<-permutest(rda(Q2~ENV4+Condition(ENV2)+Condition(ENV3)+Condition(ENV1)),permutation=1000,model="full")
		Se4<-sort(ENV4.sig$F.perm)
		if(length(Se4[Se4>ENV4.sig$F.0])==0){
			result3[,5]<-c("<0.001")
		}
		else{	
			result3[,5]<-length(Se4[Se4>ENV4.sig$F.0])/length(ENV4.sig$F.perm)		
		}
		#result3<-list(Pw,Pe1,Pe2,Pe3,Pe4)
		#names(result3)<-c("whole.sig","ENV1.sig","ENV2.sig","ENV3.sig","ENV4.sig")
	}

return(result3)
}

######################################################################
