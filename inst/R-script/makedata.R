### Preparation of the example data.
library(affy)
library(OrderedList)

load("/nfs/compdiag/molgen/gene_expression/analyses/Singh03/expression_set.rdat")
prostate <- expr.set
load("/project/gene_expression/analyses/Huang03/expression_set.rdat")
breast <- expr.set
rm(expr.set)


set.seed(123)
nn <- rownames(exprs(breast))[sample(1000)]

map <- data.frame(prostate=nn,breast=paste(nn,"_B",sep=""))

e <- exprs(prostate)[nn,]
rownames(e) <- map$prostate
p <- pData(prostate)
x <- !is.na(p$outcome)
p <- p[x,]
e <- e[,x]
pp <- as.list(names(p))
names(pp) <- names(p)
pdata <- new("phenoData", pData=p, varLabels=pp)
prostate <- new("exprSet", exprs=e, phenoData=pdata)

e <- exprs(breast)[nn,]
rownames(e) <- map$breast
p <- pData(breast)
x <- !is.na(p$Risk)
p <- p[x,]
e <- e[,x]
x1 <- rownames(p)[which(p$Risk=="high")][1:15]
x2 <- rownames(p)[which(p$Risk=="low")][1:15]
x <- c(x1,x2)
p <- p[x,]
e <- e[,x]
pp <- as.list(names(p))
names(pp) <- names(p)
pdata <- new("phenoData", pData=p, varLabels=pp)
breast <- new("exprSet", exprs=e, phenoData=pdata)

OL.data <- list(breast=breast,prostate=prostate,map=map)
save(OL.data,file="../../data/OL.data.rda",compress=T)


a <- prepareData(
                 list(data=breast,name="breast",var="Risk",out=c("high","low"),paired=FALSE),
                 list(data=prostate,name="prostate",var="outcome",out=c("Rec","NRec"),paired=FALSE),
                 mapping=map
                 )
OL.result <- OrderedList(a)

save(OL.result,file="../../data/OL.result.rda",compress=T)





