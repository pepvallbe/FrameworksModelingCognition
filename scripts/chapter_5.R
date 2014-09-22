#############################################################################################################
#TITLE: FRAMEWORKS FOR MODELING COGNITION AND DECISIONS IN INSTITUTIONAL ENVIRONMENTS: A DATA-DRIVEN APPROACH
#AUTHOR: JOAN-JOSEP VALLBÉ
#REPLICATION MATERIALS
#CHAPTER 5: REPRESENTING ORGANIZATIONAL UNCERTAINTY 
#SEPTEMBER 18, 2014
##############################################################################################################

#The document assumes you have the data stored in a subdirectory called data/
#The document assumes you have a subdirectory named graphs/
#If you install the ProjectTemplate R package it will create all the subdirectories you need



library(topicmodels)
library(tm)
library(slam)
library(FactoMineR)
library(vegan)
library(xtable)
library(SDMTools)
library(ClustOfVar)
library(TTR)
library(nnet)
library(stargazer)
library(effects)


#LOAD THE CORPUS CONTAINING 3 TYPES OF DOCUMENTS (ON-CALL, CRIMINAL, CIVIL)
load("data/corpus3.RData")

#DOCUMENT-TEXT MATRIX
dtm <- DocumentTermMatrix(corpus, control = list(weighting = weightTf))

#dtidf <- DocumentTermMatrix(corpus, control = list(weighting = weightTfIdf))

dtm$dimnames$Docs <- substr(dtm$dimnames$Docs,1,7)

dim(dtm)#we have 330 documents and 5149 words


summary(col_sums(dtm))
term_tfidf <- tapply(dtm$v/row_sums(dtm)[dtm$i], dtm$j, mean) * log2(nDocs(dtm)/col_sums(dtm > 0))

summary(term_tfidf)


dtm <- dtm[,term_tfidf>=0.035]#we put the threshold between the mean and the median tfidf, so that very frequent terms do not pop up 

dtm <- dtm[row_sums(dtm)>0,]

dim(dtm)#330 documents, 1912 terms


summary(col_sums(dtm))



########################################
#FIG 5.1: DISTRIBUTION OF TERM FREQUENCY
########################################

w <- as.matrix(dtm)
wf <- apply(w,1,sum)

summary(wf)

wf <- sort(wf,decreasing=FALSE)

pdf("graphs/hist_freq_corpora.pdf")
hist(wf,breaks=20,
     ylab="Number of documents",
     xlab="Term frequency" ,
     main="")
abline(v=mean(wf),lty=2)
dev.off()


#####################################
#HIERARCHICAL CLUSTERING
#####################################

mat <- as.matrix(dtm)

d <- vegdist(mat)

h <- hclust(d,method="ward")

names(h)

h$height

plot(h,cex=0.5)

pdf("graphs/dendrogram_all.pdf",width=15)
plot(h,cex=0.5,
     #horiz=TRUE,
     main="",
     xlab="Clusters of documents",
     sub="",
     yaxt="n")
axis(2, at=1:5,las=2)
dev.off()

colLab <- function(n) {
  if(is.leaf(n)) {
    a <- attributes(n)
#    attr(n, "label") <- substr(a$label,1,2)             #  change the node label 
    attr(n, "nodePar") <- c(a$nodePar, lab.cex = 0.2) #   change the node size
  }
  n
}

dendro <- as.dendrogram(h)
dendro <- dendrapply(dendro, colLab)

plot(dendro,horiz=TRUE)

pdf("graphs/dendrogram_all_horiz.pdf",width=15)
#par(cex=0.4,font=1)
plot(h, cex=0.3,
     main="",
     sub="",
     xlab="Clusters")
dev.off()


clust <- cutree(h,k=3)


################################
#MULTIDIMENSIONAL SCALING
################################

fitmds <- cmdscale(d,eig=TRUE,k=2)

corp <- substr(row.names(fitmds$points),7,8)

dat <- data.frame(fitmds$points,corp)

x <- cbind(dat,clust)

###################################################################
#FIG 5.3: PLOT OF THE POSITION OF ALL DOCUMENTS ON 2 MDS DIMENSIONS
###################################################################

pdf("graphs/mds_3_corpora.pdf",width=10)
plot(x$X1,x$X2,col=x$corp,pch=c(20,2,3,4)[x$corp],
     xlab="Dimension 1",
     ylab="Dimension 2")
legend("topright",col=x$corp,legend=c("civil", "on-call", "criminal"),pch=c(20,2,3)[x$corp])
dev.off()


###################################################################
#FIG 5.4: PLOT OF THE POSITION OF ALL DOCUMENTS ON 2 MDS DIMENSIONS,
#WITH ELLIPSES FROM CLUSTERING
###################################################################

pdf("graphs/mds_3_corpora_ellipse.pdf",width=10)
plot(x$X1,x$X2,col=x$corp,pch=c(20,2,3,4)[x$corp],
     xlab="Dimension 1",
     ylab="Dimension 2")
ordiellipse(fitmds,clust,label=clust)
legend("topright",col=x$corp,legend=c("civil", "on-call", "criminal"),pch=c(20,2,3)[x$corp])
dev.off()


########################################
#TABLE 5.4: CLUSTERS OF DOCUMENTS (ROW)
########################################

#Row percentages of the clusters in MDS
xtable(round(prop.table(table(x$corp,x$clust),1)*100,2))


##########################################
#TABLE 5.5: CLUSTERS OF DOCUMENTS (COLUMN)
##########################################

#Column percentages of the clusters in MDS
xtable(round(prop.table(table(x$corp,x$clust),2)*100,2))


############################################
#FIG 5.5: BOXPLOT DISTRIBUTION OF MDS SCORES
############################################

pdf("graphs/boxplotMDS.pdf")
boxplot(x$X1~x$clust, range=0,
        xlab="Clusters",
        ylab="Dimension 1")
#abline(h=0,lty=2)
dev.off()

###############################################################
#TABLE 5.6: OLS REGRESSION MODEL TESTING CORRESPONDENCE BETWEEN
#CLUSTERS IN DATA AND MDS SCORE IN FIRST DIMENSION
###############################################################

#We test whether there is a correspondence between the clusters of the data and the score in the first dimension o the MDS
mod <- lm(X1~as.factor(clust),data=x)

#We test whether there is a correspondence between the type of document (type of corpus) and the score in the first dimension of MDS
mod2 <- lm(X1~as.factor(corp),data=x)

mod2b <- update(mod2,.~. + as.factor(clust))

an2 <- anova(mod2)

#export results in ascii
stargazer(mod2,mod,mod2b,type="text")

#export results in LaTeX
stargazer(mod2,mod,mod2b,type="latex")



#########################################################
#REPLICATION WITH A CORPUS OF ONLY TWO TYPES OF DOCUMENTS
#########################################################


#LOAD THE TEXT CORPUS WITH ONLY TWO TYPES OF DOCUMENT

load("data/corpus2.RData")

#DOCUMENT-TEXT MATRIX
dtm <- DocumentTermMatrix(corpus, control = list(weighting = weightTf))

#dtidf <- DocumentTermMatrix(corpus, control = list(weighting = weightTfIdf))
dtm$dimnames$Docs <- substr(dtm$dimnames$Docs,1,7)

dim(dtm)#we have 219 documents and 4288 words


summary(col_sums(dtm))
term_tfidf <- tapply(dtm$v/row_sums(dtm)[dtm$i], dtm$j, mean) * log2(nDocs(dtm)/col_sums(dtm > 0))

summary(term_tfidf)


dtm <- dtm[,term_tfidf>=0.035]#we put the threshold between the median and the mean tfidf, so that very frequent terms do not pop up

dtm <- dtm[row_sums(dtm)>0,]

dim(dtm)#219, 1375 terms

summary(col_sums(dtm))

#dtm <- dtm[,col_sums(dtm)>1]

##########################
#HIERARCHICAL CLUSTERING
##########################

mat <- as.matrix(dtm)

d <- vegdist(mat)

h <- hclust(d,method="ward")

names(h)


clust <- cutree(h,k=2)

##########################
#MULTIDIMENSIONAL SCALING
###########################

fitmds <- cmdscale(d,eig=TRUE,k=2)

corp <- substr(row.names(fitmds$points),7,8)

dat <- data.frame(fitmds$points,corp)

x <- cbind(dat,clust)


#############################################################
#FIG 5.6: REPLICATION OF MDS WITH ONLY TWO TYPES OF DOCUMENTS
#############################################################

pdf("graphs/mds_2_corpora.pdf",width=10)
plot(x$X1,x$X2,col=x$corp,pch=c(20,3)[x$corp],
     xlab="Dimension 1",
     ylab="Dimension 2")
legend("topright",col=x$corp,legend=c("on-call", "criminal"),pch=c(20,3)[x$corp])
abline(h=0,lty=2)
#abline(v=0,lty=2)
dev.off()




####################################################
#TOPIC MODELS: LATENT DIRICHLET ALLOCATION (LDA)
####################################################

#define the text corpora with 3 types of documents

load("data/corpus3.RData")


dtm <- DocumentTermMatrix(corpus, control = list(weighting = weightTf))

#dtidf <- DocumentTermMatrix(corpus, control = list(weighting = weightTfIdf))

dtm$dimnames$Docs <- substr(dtm$dimnames$Docs,1,7)

dim(dtm)#we have 330 documents and 5149 words


summary(col_sums(dtm))
term_tfidf <- tapply(dtm$v/row_sums(dtm)[dtm$i], dtm$j, mean) * log2(nDocs(dtm)/col_sums(dtm > 0))

summary(term_tfidf)


dtm <- dtm[,term_tfidf>=0.035]#we put the threshold just above the median tfidf, so that very frequent terms do not pop up

dtm <- dtm[row_sums(dtm)>0,]

dim(dtm)#330 docs / 1912 terms

summary(col_sums(dtm))





#########################################################
#IMPLEMENTATION OF LDA TOPIC MODELS FOR the 3 corpora
#########################################################


#THIS IS THE SIMULATION OF TOPIC MODELS UP TO 100 TOPICS
#ALERT: THIS CAN TAKE SOME TIME TO COMPUTE

SEED <- 2010
n <- 100
corp.lda <- list()
for (i in 2:n) {
corp.lda[[i]] <- LDA(dtm,k=i,
    control=list(seed=SEED))}


#THIS IS TO JUST FIT THE TOPIC MODEL WITH 10 TOPICS

SEED <- 2010

corp.lda <- LDA(dtm,k=10,
     control=list(seed=SEED))


###############################################
#PLOT RELEVANT TERM ASSOCIATIONS FOR EACH TOPIC
###############################################

terms(corp.lda,10)

c <- corp.lda

t <- terms(c,10)


#PERPLEXITY
#only when corp.lda is the list containing the topic models up to 100 topics
perp <- matrix(data=NA,nrow=length(corp.lda),ncol=1)
for (i in 2:length(corp.lda)) {
    perp[i,] <- perplexity(corp.lda[[i]])
    }

########################################
#FIG 5.8 PERPLEXITY OF 100 TOPIC MODELS
########################################

pdf("graphs/perplexity_corp_100topics.pdf",width=10)
plot(perp,type="l",
     ylim=c(0,500),
     ylab="Perplexity of the model",
     xlab="Number of topics")
dev.off()


##########################################
#FIT 5.9 DIFF IN PERPLEXITY OF TOPIC MODEL
##########################################

mom <- momentum(perp)

pdf("graphs/momentum_corp.pdf",width=10)
plot(mom,type="l",
     xlab="Number of topics",
     ylab="Differences in perplexity of topic model")
abline(h=0,lty=2)
dev.off()


######################################
#FIG 5.10 RATE OF CHANGE IN PERPLEXITY
######################################


roc <- ROC(perp,type="continuous")

pdf("graphs/roc_corp.pdf",width=10)
plot(roc,type="l",
     xlab="Number of topics",
     ylab="Rate of change in perplexity of topic model",
     xaxt="n")
abline(h=0,lty=2)
abline(v=9,lty=2)
axis(1,at=c(0,9,20,40,60,80,100),labels=c("0","9","20","40","60","80","100"))
dev.off()


##########################################
#TABLE 5.8 DISTRIBUTION OF DOCS PER TOPIC
##########################################

c <- corp.lda

top <- as.data.frame(table(topics(c)))

perc <- round((top$Freq/sum(top$Freq))*100,2)

top <- data.frame(top,perc)

xtable((top[,-1]))


############################################
#FIG 5.11 DISTRIBUTION OF TOPICS IN 20 DOCS
############################################

post <- posterior(c)$topics

pdf("graphs/distribution_topics_documents.pdf",width=10)
par(mfrow=c(4,5))
for (i in 1:20){
     barplot(post[i,],
             ylim=c(0,1),
             main= i ,
             xlab="Topics",
             ylab="Probability")
    }
dev.off()


##################################################################
#TABLE 5.9: LIST OF THE MOST RELEVANT TERMS IN EACH TOPIC, WITH P
##################################################################

d <- (posterior(c)$terms)

top1 <- sort(d[1,],decreasing=TRUE)[1:20]
top2 <- sort(d[2,],decreasing=TRUE)[1:20]
top3 <- sort(d[3,],decreasing=TRUE)[1:20]
top4 <- sort(d[4,],decreasing=TRUE)[1:20]
top5 <- sort(d[5,],decreasing=TRUE)[1:20]
top6 <- sort(d[6,],decreasing=TRUE)[1:20]
top7 <- sort(d[7,],decreasing=TRUE)[1:20]
top8 <- sort(d[8,],decreasing=TRUE)[1:20]
top9 <- sort(d[9,],decreasing=TRUE)[1:20]
top10 <- sort(d[10,],decreasing=TRUE)[1:20]

tops <- data.frame(names(top1),top1,
                   names(top2),top2,
                   names(top3),top3,
                   names(top4),top4,
                   names(top5),top5,
                   names(top6),top6,
                   names(top7),top7,
                   names(top8),top8,
                   names(top9),top9,
                   names(top10),top10,row.names=NULL)

xtable(tops[,11:20],digits=2)


##########################################
#FIG. 5.12-5.20 TERM CORRELATION NETWORKS
##########################################

pdf("graphs/corr_terms_T1.pdf")
plot(dtm, terms(corp.lda,15)[,1],corThreshold=0.25)
dev.off()

pdf("graphs/corr_terms_T2.pdf")
plot(dtm, terms(corp.lda,15)[,2],corThreshold=0.2)
dev.off()

pdf("graphs/corr_terms_T3.pdf")
plot(dtm, terms(corp.lda,15)[,3],corThreshold=0.25)
dev.off()

pdf("graphs/corr_terms_T4.pdf")
plot(dtm, terms(corp.lda,15)[,4],corThreshold=0.2)
dev.off()

pdf("graphs/corr_terms_T5.pdf")
plot(dtm, terms(corp.lda,15)[,5],corThreshold=0.25)
dev.off()

pdf("graphs/corr_terms_T6.pdf")
plot(dtm, terms(corp.lda,15)[,6],corThreshold=0.25)
dev.off()

pdf("graphs/corr_terms_T7.pdf")
plot(dtm, terms(corp.lda,15)[,7],corThreshold=0.2)
dev.off()

pdf("graphs/corr_terms_T8.pdf")
plot(dtm, terms(corp.lda,15)[,8],corThreshold=0.2)
dev.off()

pdf("graphs/corr_terms_T9.pdf")
plot(dtm, terms(corp.lda,15)[,9],corThreshold=0.25)
dev.off()

pdf("graphs/corr_terms_T10.pdf")
plot(dtm, terms(corp.lda,15)[,10],corThreshold=0.2)
dev.off()



#Distribution of topics per one document

post <- posterior(c)$topics

p <- PCA(as.matrix(post), graph=FALSE)

################################################################
#TABLE 5.11 CORRELATION OF EACH TOPIC TO THE FIRST PCA COMPONENT
################################################################

dimdesc(p)[[1]]


################################################################
#FIG 5.21 CLASSIFICATION OF TOPICS ALONG THE FIRST PCA COMPONENT
################################################################

pdf("graphs/pca_topics.pdf")
plot(p$var$coord[,1:2],
     type="p",
     pch=20,
     xlim=c(-.9,.9),
     ylim=c(-.7,1.2),
     yaxt="n",
     xlab="First principal component",
     ylab="")
abline(v=0,lty=2)
#points(p$ind$coord[,1:2], pch=20,col="grey")
text(p$var$coord[,1:2],as.character(rownames(p$var$coord)),cex=1.1,pos=4)
arrows(0.2,.9, .8, .9, length = 0.15, angle = 10,code=2,
       lty=1,lwd=par("lwd"),col="azure4")
text(0.5,1,"behavioral",col="azure4")
arrows(-0.2,.9, -.8, .9, length = 0.15, angle = 10,code=2,
       lty=1,lwd=par("lwd"),col="azure4")
text(-0.5,1,"theoretical",col="azure4")
dev.off()



#####################################
#TOPIC CLASSIFICATION
#####################################


corp <- substr(row.names(post),7,8)

dat <- data.frame(post,corp)

maxims <- apply(dat[,1:10],1,max)

#extract the variable with a highest value
dat$tops <- apply(dat[, -11], 1, function(x) which(x == max(x)))
dat$weight <-  apply(dat[,1:10],1,max)


#Our hypothesis is that Topics 1, 3, 7, 9 are technical, and 2, 5 ,8 are behavioral. Topic 6 is neutral. We recode the top variable: behavioral=1, neutral=0, technical=-1  

dat$thtop[as.factor(dat$tops)=="1" |
          as.factor(dat$tops)=="3" |
          as.factor(dat$tops)=="7" |
          as.factor(dat$tops)=="9" |
          as.factor(dat$tops)=="10"] <- 1
dat$thtop[as.factor(dat$tops)=="2" |
          as.factor(dat$tops)=="5" |
          as.factor(dat$tops)=="8" ] <- -1
dat$thtop[as.factor(dat$tops)=="4" |
          as.factor(dat$tops)=="6"] <- 0

#Classification suggested by PCA results


dat$pcatop[as.factor(dat$tops)=="1" |
          as.factor(dat$tops)=="3" |
          as.factor(dat$tops)=="7" |
          as.factor(dat$tops)=="9"] <- 1
dat$pcatop[as.factor(dat$tops)=="2" |
          as.factor(dat$tops)=="8" |          
          as.factor(dat$tops)=="5" ] <- -1
dat$pcatop[as.factor(dat$tops)=="4" |
          as.factor(dat$tops)=="6" |
          as.factor(dat$tops)=="10"] <- 0


##################################################################
#TABLE 5.12 CROSSTABS BETWEEN THEORETICAL AND DATA CLASSIFICATIONS
##################################################################

xtabs(~thtop + pcatop,data=dat)

########################################################
#TABLE 5.13 CROSSTABS TYPE OF CORPUS WITH CLASSIFICATION
########################################################

#theoretical model
round(prop.table(table(dat$corp,dat$thtop),1)*100,2)

#data driven model
round(prop.table(table(dat$corp,dat$pcatop),1)*100,2)



#########################################################
#TABLE 5.14 MULTINOMIAL REGRESSION MODELS
#########################################################

#model1: theoretical classification

mod1 <- multinom(thtop~corp,  data=dat)

#model2: data-driven classification

mod2 <- multinom(pcatop~corp,  data=dat)

#Export results ascii
stargazer(mod1,mod2,type="text",digits=2)

#Export results to LaTeX
stargazer(mod1,mod2,type="latex",digits=2)


################################################
#TABLE 5.15 PREDICTED PROBABILITIES
################################################

dcorpus <- data.frame(corp = c("c", "p", "g"))

#Theoretical model
round(predict(mod1, newdata = dcorpus, "probs"),2)

#Data-driven model
round(predict(mod2, newdata = dcorpus, "probs"),2)


#############################################
#FIG 5.22 PLOT OF PREDICTED PROBABILITIES
#############################################

pdf("graphs/effects_theoretical_model.pdf")
plot(allEffects(mod1),
     xlab="Type of Corpus",
     ylab="Probability",
     main="Theoretical model")
dev.off()

pdf("graphs/effects_data_model.pdf")
plot(allEffects(mod2),
     xlab="Type of Corpus",
     ylab="Probability",
     main="Data-driven model")
dev.off()
