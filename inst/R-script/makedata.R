### Preparation of the example data.
library(affy)
library(OrderedList)

source("/nfs/compdiag/user/koc16861/bioconductor/OrderedList/R/prepareData.R")
source("/nfs/compdiag/user/koc16861/bioconductor/OrderedList/R/OrderedList.R")

load("/nfs/compdiag/molgen/gene_expression/analyses/Singh03/expression_set.rdat") 
cat("class(expr.set_prostate) ::",class(expr.set),"\n") # AffyBatch
expr.set.prostate<-expr.set
cat("class(expr.set.prostate@phenoData)::", class(expr.set.prostate@phenoData),"\n")
prostate <- updateObject(expr.set)
cat("class(prostate)[update] ::",class(prostate),"\n") # AffyBatch

load("/nfs/compdiag/molgen/gene_expression/analyses/Huang03/expression_set.rdat")
cat("class(expr.set_breast) ::",class(expr.set),"\n") # AffyBatch
expr.set.breast<-expr.set
cat("class(expr.set.breast@phenoData)::", class(expr.set.breast@phenoData),"\n")
breast <- updateObject(expr.set)
cat("class(breast)[update] ::",class(breast),"\n") # AffyBatch

rm(expr.set)


set.seed(123)
nn <- rownames(intensity(breast))[sample(1000)]

map <- data.frame(prostate=nn,breast=paste(nn,"_B",sep=""))

e <- intensity(prostate)[nn,] # class(prostate):= AffyBatch

rownames(e) <- map$prostate
p <- pData(prostate)	# phenoData (rows: samples ; col: covariates)
x <- !is.na(p$outcome)	# extract content of field 'outcome' of phenoData.


p <- p[x,] #take only those samples, where !is.NA is TRUE(=any value is available), e.g. p[1,1]=NA => Sample 1 isSkipped
e <- e[,x] 	# skip the same samples for 'prostate'-intensities


pp <- data.frame(labelDescription=names(p)) # pp:=metadata (rows: covariates, 1 column: labelDescription)
rownames(pp) <- colnames(p)
pdata <- new("AnnotatedDataFrame", data=p, varMetadata=pp)
prostate <- new("ExpressionSet", exprs=e, phenoData=pdata)

e <- intensity(breast)[nn,]

rownames(e) <- map$breast
p <- pData(breast)
x <- !is.na(p$Risk)
p <- p[x,]
e <- e[,x]

x1 <- rownames(p)[which(p$Risk=="high")][1:15] #which indices are TRUE? | taken from breast data
x2 <- rownames(p)[which(p$Risk=="low")][1:15]

x <- c(x1,x2)
p <- p[x,]
e <- e[,x]

pp <- data.frame(labelDescription=names(p))
rownames(pp) <- colnames(p)
pdata <- new("AnnotatedDataFrame", data=p, varMetadata=pp)
breast <- new("ExpressionSet", exprs=e, phenoData=pdata)


OL.data <- list(breast=breast,prostate=prostate,map=map)
save(OL.data,file="../../data/OL.data.rda",compress=T)

a <- prepareData(		                 		
		list(data=breast,name="breast",var="Risk",out=c("high","low"),paired=FALSE),
		list(data=prostate,name="prostate",var="outcome",out=c("Rec","NRec"),paired=FALSE),	             			mapping=map
               )

#eset1<-list(data=breast,name="breast",var="Risk",out=c("high","low"),paired=FALSE)

OL.result <- OrderedList(a)
save(OL.result,file="../../data/OL.result.rda",compress=T)





