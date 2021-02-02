# This is my example with notes 


file <- "TCGA_breast_cancer_LumA_vs_Basal_PAM50.tsv" #variable type: string
first10 <- c('NAT1','BIRC5','BAG1','BCL2','BLVRA','CCNB1','CCNE1','CDC6','CDC20','CDH3')
nfold <- 3
header <- scan(file, nlines = 1, sep="\t", what = character()) #use the function of "scan" to read the file, and read only one line
#print (header)
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
names(data) <- header
#print (data)
#dim(data) #dimension of the data--->rows,columns
#data$sample_id #show column labeled “sample_id” – gene names
#data$sample_id %in% c("BCL2","NAT1") #apply operater %in% to filter 2 genes, c()works for make the genes a list [c:concatenate]
#as.character(data$sample_id) #a better way to see genes
#data$sample_id == "BCL2" #filer the sample id to be BCL2
#data[data$sample_id=="BCL2",] #All samples from gene "BCL2". Don't forget the comma
#data[,c(6)] #data of column 6
#data[c(6)] #data of colume 6, easier to watch

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character()) #read the line of "class"
#data[,header2=='Luminal A'] #filer on columns

LumA <- data[data$sample_id %in% first10,header2=='Luminal A']
Basal <- data[data$sample_id %in% first10,header2=='Basal-like']

LumA_groups <- split(colnames(LumA), sample(1:nfold, ncol(LumA), replace=T))
Basal_groups <- split(colnames(Basal), sample(1:nfold, ncol(Basal), replace=T))
  
result <- array()
  
for (test_group in 1:nfold) {
    
    testA <- LumA[,colnames(LumA) %in% unlist(LumA_groups[test_group])]
    testB <- Basal[,colnames(Basal) %in% unlist(Basal_groups[test_group])]
    
    trainingA <- LumA[,!(colnames(LumA) %in% unlist(LumA_groups[test_group]))]
    trainingB <- Basal[,!(colnames(Basal) %in% unlist(Basal_groups[test_group]))]
    
    centroidA <- rowMeans(trainingA)
    centroidB <- rowMeans(trainingB)
    
    misclassifiedA <- sum(sapply(testA, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))>0 }))
    misclassifiedB <- sum(sapply(testB, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))<0 }))
    
    result[test_group] <- (misclassifiedA+misclassifiedB)/(ncol(testA)+ncol(testB))
}
  
mean(result)
sd(result)

