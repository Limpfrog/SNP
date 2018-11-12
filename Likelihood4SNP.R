
# Ali Amiryousefi
#
# ali.amiryousefi@helsinki.fi
#
# Maximum likelihood score for SNP data
#
# The relative existence probability/score for a species.



prob.Filler<- function(fasta, prob=TRUE){
	nr<- length(fasta)
	nc<- length(fasta[[1]])
	p<- matrix(0, nr, nc)
	n<- matrix(0, nr, nc)
	for (i in 1:nr){
		for (j in 1:nc){
			n[i, j]<- fas[[i]][j]
		}
	}
	for (i in 1:nc){
		sumn<- sum(n[,i]=="a")+sum(n[,i]=="c")+sum(n[,i]=="g")+sum(n[,i]=="t")
		probA<- sum(n[,i]=="a")/sumn
		probC<- sum(n[,i]=="c")/sumn
		probG<- sum(n[,i]=="g")/sumn
		probT<- sum(n[,i]=="t")/sumn 
	col.prob<- c(probA, probC, probG, probT)
	for (j in 1:nr){
		if (n[j, i]=="a"){p[j,i]<- probA}
		else if (n[j, i]=="c"){p[j,i]<- probC}
		else if (n[j, i]=="g"){p[j,i]<- probG}
		else if (n[j, i]=="t"){p[j,i]<- probT}
		else {p[j,i]<- NA}
		}	
	}
	zero<- !logical(nc)		#removing the columns with only 1 outputs. Monocharacter columns!
	for (j in 1:nc){
		if (sum(p[,j]==1, na.rm=TRUE)>0){
			zero[j]<- FALSE
		}
	}
	p<- subset(p, select=zero)
	for (i in 1:nrow(p)){# point estimating the missing values with the 1/4 probablity of that state being any nucleotide (note this is different than GC content since the overal probablity of each state was apporved to be close to 25%)
		for (j in 1:ncol(p)){
			if (is.na(p[i,j])){p[i,j]<-0.25}
		}
	}
	likelihood<- numeric(length(p[,1]))
	mll<- numeric(length(p[,1])) #minus ln likelihood
	mll10<- numeric(length(p[,1])) #minus log_10 likelihood
	for (i in 1:length(p[,1])){
		mll[i]<- sum(-log(p[i,]))
		mll10[i]<- sum(-log10(p[i,]))
		}
	out<- data.frame("speecies"=names(fasta), "nMLL"= mll, "MLL10"=mll10)
	return(out)
}# returning the table of species names and the two columns for (minus) natural and log likelihood of each sequence as the bigger the value the smaller the probability of that species hence more exotic in comparisons to the others. 

#read the data in fas
#fas<- read.fasta("~/Desktop/SNP_convert_Final.fasta")


p<- prob.Filler(fas)