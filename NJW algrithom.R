#NJW Algorithm 
#toy data simulation
#load library for multivariate normal random data
library(MASS)

#specify the means
mu1 <- c(35,10)
mu4 <- c(10,55)

#specify the covariance matrix (multivariate version of variance)
sigma.matrix <- matrix(c(1,0,0,1),2,2,byrow = FALSE)

#check if the covariance matrix is positive definite
#we want all of the eigenvalues to be positive
eigen(sigma.matrix)

#make the data randomly from 2 populations
set.seed(1000)

gaussian1 <- mvrnorm(n = 200, mu=mu1, Sigma = sigma.matrix)
gaussian4 <- mvrnorm(n = 200, mu=mu4, Sigma = sigma.matrix)

my.gaussian.data <- rbind(gaussian1,gaussian4) #400*2
head(my.gaussian.data)
plot(my.gaussian.data)
#go to line 124 for scalable NJW

#####Original NJW#####

#firstly apply L2 normalization for A(cos similiarty)

L2.Normalization=function(x){
  matrix.normalization=normalize.rows(x)
  return(matrix.normalization)
}



#step1: compute W matrix W=A %*% t(A)
W.matrix<-function(x){
 w=tcrossprod(x)
 w[which(w<0)]=0
 return(w)
}


#step2: compute D
D.matrix.original<-function(x,alpha){
  n1=nrow(x)
  D=rowSums(x)
  labels=rep(T,n1)
  if(alpha>0){
    D.keep=D[-order(D)[1:round(n1*alpha)]]
    D.out=D.keep
    A.out=x[-order(D)[1:round(n1*alpha)],]
    labels[order(D)[1:round(n1*alpha)]]=F
  }else{
    A.out=x
    D.out=D}
  return(list(D=D.out,A=A.out,label=labels))
}


#setp3: compute W_tilda
W_tilda<-function(W,D){
  D.sqrt=as.numeric(D**(-0.5))
  w_tilda<-Diagonal(x=D.sqrt) %*% W %*% Diagonal(x=D.sqrt) 
  return(w_tilda)
}


#step4:eigen decomposition to w_tilda compute U
U_matrix<-function(x,k){
  u=eigs(x,k,which = "LM")$vectors
  return(u)
}


#step5: Uk row normalizations
#L2 normalization
k=2
new.U=U[,1:k]
v=L2.normalization(new.U)

#original NJW function
original.NJW<-function(A, alpha,k){
  A.nor<-L2.Normalization(A)
  W=W.matrix(A.nor)
  D=D.matrix.original(W,alpha)$D
  W_tilda=W_tilda(W,D)
  U=U_matrix(W_tilda,k)
  V=L2.Normalization(U)
  cluster<- kmeans(x = V, centers = k, iter.max = 10, nstart = 10)$cluster
  return(cluster)
}

t1=proc.time()
a<-original.NJW(usps_data,0,10)
proc.time()-t1
accuracy(usps_labels,a)
###Scalable NJW########

#step0: normalize A

library(Matrix)
library(wordspace)


#L2.Normalization=function(x){
#  row.norm=apply(x,1,norm,type="2")
#  matrix.normalization=x/row.norm
#  return(matrix.normalization)
#}

x=usps_data

L2.Normalization=function(x){
  matrix.normalization=normalize.rows(x)
  return(matrix.normalization)
}

q<-L2.Normalization(x)
#step1: calculate D
D.matrix=function(x,alpha){
  n1=nrow(x)
  one=matrix(1,nrow = n1)
  D=x%*%(t(x)%*%one)-one
  labels=rep(T,n1)
  if(alpha>0){
    D.keep=D[-order(D)[1:round(n1*alpha)]]
    D.out=D.keep
    A.out=x[-order(D)[1:round(n1*alpha)],]
    labels[order(D)[1:round(n1*alpha)]]=F
  }else{
    A.out=x
    D.out=D}
  return(list(D=D.out,A=A.out,label=labels))
}

D<-D.matrix(q,0.01)$D
plot(sort(1/D,decreasing=T),col=ifelse(a>0.00028,'red','blue'),xlab='',ylab='',xaxs='i')
axis(1,at=seq(0,9000,by=1000))



a=sort(1/D,decreasing=T)
summary(a)
A<-D.matrix(q,0)$A
a<-D.matrix(usps_data,0)$D
a<-a^(-1)
plot(sort(a))


#Step 2: calculate A_tilda
A.tilda=function(D,A){
 D.sqrt=as.numeric(D**-0.5)
 A.tilda=Diagonal(x=D.sqrt) %*% A
 return(A.tilda)
}
#a<-A.tilda(D,A)


#install.packages("RSpectra")
library(RSpectra)
#step 3: Apply SVD on A
A.svd<-function(A,k){
  out=svds(A,k)
  return(out$u)
}


#accurarcy: from group2
library(clue)
accuracy <- function(truelabels, clusters) {
  #Hungarian Algorithm
  
  # labels from cluster A will be matched on the labels from cluster B
  
  minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
    require(clue)
    idsA <- unique(clusteringA)  # distinct cluster ids in a
    idsB <- unique(clusteringB)  # distinct cluster ids in b
    nA <- length(clusteringA)  # number of instances in a
    nB <- length(clusteringB)  # number of instances in b
    if (length(idsA) != length(idsB) || nA != nB) {
      stop("number of clusters do not match")
    }
    
    nC <- length(idsA)
    tupel <- c(1:nA)
    
    # computing the assignment matrix
    assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
    for (i in 1:nC) {
      tupelClusterI <- tupel[clusteringA == i]
      solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
        nA_I <- length(tupelA_I)  # number of elements in cluster I
        tupelB_I <- tupel[clusterIDsB == i]
        nB_I <- length(tupelB_I)
        nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
        return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
      }, clusteringB, tupelClusterI)
      assignmentMatrix[i, ] <- solRowI
    }
    
    # optimization
    result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
    attr(result, "assignmentMatrix") <- assignmentMatrix
    return(result)
  }
  
  test <- minWeightBipartiteMatching(clusters, truelabels)
  
  predicted = NULL; 
  predicted <- rep(NA, length(clusters))
  for(i in 1:length(test)) {
    predicted[which(clusters == i)] <- test[i] 
  }
  
  table <- table(predicted, truelabels)
  accuracy <- (sum(diag(table))/length(truelabels))
  classaccuracy <- vector()
  colsums <- colSums(table)
  for(i in 1:length(test)) {
    classaccuracy[i] <- table[i,i]/colsums[i]
  }
  
  return(list("accuracy" = accuracy, "classaccuracy" = classaccuracy, "table" = table,
              "mapping" = test, "mappedlabels" = predicted))
}


#combine everything together
sca.NJW<-function(A,alpha,k){
  A.nor=L2.Normalization(A)
  D=D.matrix(A.nor,alpha)$D
  A=D.matrix(A.nor,alpha)$A
  Atilda=A.tilda(D,A)
  U=A.svd(Atilda,k)
  v=L2.Normalization(U)
  set.seed(1000)
  Kmeans=kmeans(x = v, centers = k, iter.max = 100, nstart = 10)
  Cluster=Kmeans$cluster

  return(Cluster)
}

#obtaining outlier label
outlier.function<-function(data,alpha){
  data.nor=L2.Normalization(data)
  label=D.matrix(data.nor,alpha)$label 
  return(label)
}

#check results on toy data
predicted_cluster=sca.NJW(my.gaussian.data,0,2)
true_label=c(rep(1,200),rep(2,200))
accuracy(predicted_cluster,true_label)
plot(my.gaussian.data,col=true_label+1,pch=true_label+1)

#simulation of more data
gaussian.data<-function(n){
new.gaussian1<-rep(gaussian1,n)
new.gaussian4<-rep(gaussian4,n)
my.gaussian.data <- cbind(new.gaussian1,new.gaussian4)
return(my.gaussian.data)
} #400*2

#head(gaussian.data(2))
#dim(gaussian.data(2))


i=1
data.list<-list()
while(i<=10){
 data.list[[i]]<-gaussian.data(i)
  i=i+1
}
str(data.list)

running.time<-function(x){
  elapse.time=NULL
for(i in 1:x){
  elapse.time[i]=system.time(sca.NJW(data.list[[i]],0,2))[2]
}
  return(elapse.time)
}

#try with 400-4000
time.vector=running.time(10)
n=400*c(1:10)
plot(n,time.vector,type="o",xlab="sample size",ylab="time",
     main="Running time")



#usps data without remove outliers
#install.packages("R.matlab")
setwd("~/Dropbox/sjsu/SPRING 2018/MATH 203/data")
library(R.matlab)
usps_mat <- readMat("usps.mat")
usps_data <- usps_mat$fea
dim(usps_data)
usps_labels <- usps_mat$gnd
length(usps_labels)
table(usps_labels)


#without removing outliers
predicted_cluster.original<-sca.NJW(usps_data,0,10)
accuracy(usps_labels,predicted_cluster.original)
#69.26%

#usps with 10% outlier removed

data.nor=L2.Normalization(usps_data)
data.outlier.removed=D.matrix(data.nor,0.02)$A
label=D.matrix(data.nor,0.02)$label


#obtaining outlier label
outlier.function<-function(data,alpha){
  data.nor=L2.Normalization(data)
  label=D.matrix(data.nor,alpha)$label 
  return(label)
}

outlier.label<-outlier.function(usps_data,0.02)
table(outlier.label)
predicted_clusters.10per=sca.NJW(usps_data,0.02,10)
accuracy(usps_labels[outlier.label],predicted_clusters.10per)
#accuracy=0.7247

a=rep(0,10)
for(i in 1:11){
  outlier.percent<-seq(0,0.1,by=0.01)
  data.nor=L2.Normalization(usps_data)
  label=D.matrix(data.nor,outlier.percent[i])$label
  predicted_clusters=sca.NJW(usps_data,outlier.percent[i],10)
  a[i]<-accuracy(usps_labels[label],predicted_clusters)$accuracy
}
a.vector<-unlist(a,use.names = FALSE)
outlier.percent<-seq(0,0.1,by=0.01)
plot(outlier.percent,a.vector,type='o',ylab = 'Accuracy',
     main = "Outlier Removal with USPS Data",yaxt='n')
axis(2,at=pretty(a.vector),lab=pretty(a.vector)*100,las=TRUE)


#KNN with outliers
#install.packages('class')
library(class)

#try 2% outliers first
data.nor=L2.Normalization(usps_data)
dim(data.nor)
fitted.data<-D.matrix(data.nor,0.02)$A #train
dim(fitted.data)
outlier.label<-outlier.function(usps_data,0.02)
length(outlier.label)
test<-data.nor[which(outlier.label==FALSE),]
dim(test)

predicted_clusters.10per=sca.NJW(usps_data,0.02,10)

#method 1
knn.cluster<-knn(fitted.data, test, predicted_clusters.10per, k = 25, prob=F)
table(knn.cluster)
length(knn.cluster)

#Check accuracy again

y<-c(predicted_clusters.10per,knn.cluster)
usps_labels[outlier.label]
usps_labels[which(outlier.label==FALSE)]
y.true<-c(usps_labels[outlier.label],usps_labels[which(outlier.label==FALSE)])

#length(usps_labels[which(outlier.label==FALSE)])
accuracy(y.true,y)$accuracy



#run this with k=1-30
a=rep(0,35)
for(i in 1:35){
  knn.cluster<-knn(fitted.data, test, predicted_clusters.10per, k = i, prob=F)
  y<-c(predicted_clusters.10per,knn.cluster)
  y.true<-c(usps_labels[outlier.label],usps_labels[which(outlier.label==FALSE)])
  a[i]=accuracy(y.true,y)$accuracy}
a
plot(seq(1:35),a,xlab='Number of neighbors',ylab='Accuracy',yaxt='n',
     main='Different k with 2% Outlier')
axis(2,at=pretty(a),lab=pretty(a)*100,las=TRUE)
abline(v=33,col='red')

a
max(a)

#experiment with document datasets
#reuters dataset#####
reuters_mat <- readMat("reuters.mat")
reuters_data <- reuters_mat$A
#library(Matrix)
#reuters_data<-Matrix(reuters_data,sparse = TRUE)
dim(reuters_data)
reuters_labels <- reuters_mat$y
length(reuters_labels)
table(reuters_labels)
#number of class=30


r.label<-ifelse(reuters_labels<=30,TRUE,FALSE)
table(r.label)
reuters_data.new<-reuters_data[which(r.label==TRUE),]
dim(reuters_data.new) #8067
reuters_labels.new<-reuters_labels[r.label]
length(reuters_labels.new) #8067
table(reuters_labels.new)
#without removing outliers

t1=proc.time()
cluster.reuters<-sca.NJW(reuters_data.new,0,30)
length(cluster.reuters)
proc.time()-t1
table(cluster.reuters)
accuracy(reuters_labels.new,cluster.reuters)
length(reuters_labels.new)

#loop with different number of clusters
reuter.accuracy<-rep(0,30)
for(i in 1:30){
  r.label<-ifelse(reuters_labels<=i,TRUE,FALSE)
  reuters_data.new<-reuters_data[which(r.label==TRUE),]
  reuters_labels.new<-reuters_labels[r.label]
  cluster.reuters<-sca.NJW(reuters_data.new,0,i)
  reuter.accuracy[i]<-accuracy(reuters_labels.new,cluster.reuters)$accuracy
}

plot(seq(1:30),y=reuter.accuracy,type = "o",main='Accuracy for Different Clusters',
     ylab='Accuracy(Reuters)',xlab="Number of Clusters",yaxt='n')
axis(2,at=pretty(a),lab=pretty(a)*100,las=TRUE)

save(reuter.accuracy,file = 'reuter with 30 clusters.rdata')

#looping with different percent of outliers

#without KNN
reuters.outlier=rep(0,10)
for(i in 1:11){
  outlier.percent<-seq(0,0.1,by=0.01)
  data.nor=L2.Normalization(reuters_data.new)
  label=D.matrix(data.nor,outlier.percent[i])$label
  predicted_clusters=sca.NJW(reuters_data.new,outlier.percent[i],30)
  reuters.outlier[i]<-accuracy(reuters_labels.new[label],predicted_clusters)$accuracy
}
reuters.outlier
save(reuters.outlier,file = 'reuters.outliers.without KNN.rdata')

#with KNN
reuters.KNN.accuracy=rep(0,10)
library(class)
for(i in 1:11){
outlier.percent<-seq(0,0.1,by=0.01)
data.nor=L2.Normalization(reuters_data.new)
fitted.data<-D.matrix(data.nor,outlier.percent[i])$A #train
outlier.label<-outlier.function(reuters_data.new,outlier.percent[i])
test<-data.nor[which(outlier.label==FALSE),]
predicted_clusters<-sca.NJW(reuters_data.new,outlier.percent[i],30)
knn.cluster<-knn(fitted.data, test, predicted_clusters, k = 10, prob=F)
y<-c(predicted_clusters,knn.cluster)
y.true<-c(reuters_labels.new[outlier.label],reuters_labels.new[which(outlier.label==FALSE)])
reuters.KNN.accuracy[i]<-accuracy(y.true,y)$accuracy
}

length(predicted_clusters)
dim(fitted.data)

#TDT2 dataset####
TDT2_mat <- readMat("TDT2.mat")
TDT2_data <- TDT2_mat$A
#library(Matrix)
TDT2_data<-Matrix(TDT2_data,sparse = TRUE)
dim(TDT2_data)
TDT2_labels <- TDT2_mat$y
length(TDT2_labels)
table(TDT2_labels)


cluster.TDT2<-sca.NJW(TDT2_data,0,30)
table(cluster.TDT2)
accuracy(TDT2_labels,cluster.TDT2)

#plot
library(ggplot2)




#20newsgroup

newgroup_mat<-readMat("news20.mat")
newsgroup_data <- newgroup_mat$A
dim(newsgroup_data)
newsgroup_labels <- newgroup_mat$y
length(newsgroup_labels)
table(newsgroup_labels)


#single run
set.seed(1000)
newsgroup.cluster<-sca.NJW(newsgroup_data,0,20)
accuracy(newsgroup_labels,newsgroup.cluster)$accuracy
accuracy(newsgroup_labels[outlier.function(newsgroup_data,0.01)],newsgroup.cluster)

#with loops for different alpha
newsgroup.accuracy=NULL
for(i in 1:10){
alpha=seq(0.01,0.1,by=0.01) 
newsgroup.cluster<-sca.NJW(newsgroup_data,alpha[i],20)
newsgroup.accuracy[i]=accuracy(newsgroup_labels[outlier.function(newsgroup_data,alpha[i])],newsgroup.cluster)$accuracy
}
alpha=seq(0.01,0.1,by=0.01) 


plot(alpha,newsgroup.accuracy,type = 'b',xaxt='n',yaxt='n',ylab='Accuracy',xlab='Outlier Percentage')
axis(2,at=pretty(unlist(newsgroup.accuracy)),lab=pretty(unlist(newsgroup.accuracy))*100,las=TRUE)
axis(1,at=pretty(alpha),lab=pretty(alpha)*100)

#with KNN too slow so just give up on this method
library(class)
news.data<-L2.Normalization(newsgroup_data)

train<-news.data[outlier.function(newsgroup_data,alpha[i]),]
test<-news.data[-outlier.function(newsgroup_data,alpha[i]),]
clusters<-sca.NJW(news.data,alpha[i],20)
knn.cluster<-knn(train=train,test=test,cl=clusters,k=3,prob=FALSE)

#KNN with cluster 2nd try give up with too slow KNN function

knn.accuracy=NULL
alpha=seq(0.01,0.1,by=0.01) 
center=Matrix(0,nrow=20,ncol=55570)
for(j in 1:10){

newsgroup.cluster<-sca.NJW(newsgroup_data,alpha[j],20)
outlier.label<-outlier.function(newsgroup_data,alpha[j])
newsgroup.data.withoutlier<-newsgroup_data[outlier.label,]

for(i in 1:20){
center[i,]<-apply(newsgroup.data.withoutlier[which(newsgroup.cluster==i),],2,mean)}
test<-newsgroup_data[which(outlier.function(newsgroup_data,alpha[j])==FALSE),]
knn.cluster<-knn(train=center,test=test,cl=seq(1,20),k=1,prob=FALSE)
total.cluster<-c(newsgroup.cluster,knn.cluster)
true.labels<-c(newsgroup_labels[which(outlier.label==TRUE)],
               newsgroup_labels[which(outlier.label==FALSE)])
knn.accuracy[j]<-accuracy(true.labels,total.cluster)$accuracy
}
knn.accuracy
save(knn.accuracy,file='knn.accuracy.rdata')


plot(alpha,newsgroup.accuracy,type = 'b',xaxt='n',yaxt='n',col='red',pch=5,
     ylab='Accuracy',xlab='Outlier Percentage',main='Newsgroup Accuracy with Outliers')
axis(2,at=pretty(unlist(newsgroup.accuracy)),lab=pretty(unlist(newsgroup.accuracy))*100,las=TRUE)
axis(1,at=pretty(alpha),lab=pretty(alpha)*100)
lines(alpha,knn.accuracy,type='b',col='blue',pch=8)
legend('topleft',c('Without KNN','With KNN'),cex=0.5,col=c('red','blue'),lty=c(1,1))


#Mnist
Mnist_mat <- readMat("mnist.mat")
Mnist.data<-Mnist_mat$fea
dim(Mnist_mat$fea)



#pendigits
library(rrr)
data("pendigits")
dim(pendigits)
pend_mat <- readMat("pend.mat")
pend.data<-pend_mat$fea
dim(pend.data)

library(mlbench)
a<-Mnist.data[which(Mnist.data==0)]
dim(Mnist.data)
length(a)/(70000*784)

mnist.data<-Matrix(Mnist.data,sparse = T)
class(mnist.data)
