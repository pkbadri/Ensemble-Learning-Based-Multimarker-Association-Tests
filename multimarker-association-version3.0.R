## Replaces random forests with light gbm
library(lightgbm)
library(foreach)
library(doParallel)
registerDoParallel(cores=detectCores())
library(e1071)

args      =   commandArgs(trailingOnly = TRUE)
nsnps     =   as.integer(args[2])
features  =   as.integer(args[4])
samples   =   as.integer(args[6])
trait     =   as.integer(args[8])
covs      =   as.integer(args[10])
in.snp    =   as.character(args[12])
in.pheno  =   as.character(args[14])
in.covs   =   as.character(args[16])

sfold     =   as.integer(0.1*samples)
strain    =   9*sfold;
Z         =   matrix(scan(in.snp),  nrow = samples,  ncol = nsnps, byrow=T);
P         =   scan(in.pheno);
if(covs > 0){
COV       =   matrix(scan(in.covs), nrow = samples,  ncol = covs,  byrow=T);
}
reorder   =   sample(c(1:samples))





Ptemp = P; Ztemp = Z; 
if(covs > 0){COVtemp = COV;}
for(i in 1:samples){
P[i]    = Ptemp[reorder[i]];
Z[i,]   = Ztemp[reorder[i],]; 
if(covs > 0){COV[i,] = COVtemp[reorder[i],];}
}
remove(Ptemp); remove(Ztemp); 
if(covs > 0){remove(COVtemp);}





foreach (k=1:10) %dopar% {
if(k <  10)          {testsize = sfold;    trainsize = samples - testsize;}
if(k == 10)          {trainsize = strain;  testsize = samples - trainsize;}





Ztrain    =  matrix(rep(0, (nsnps*trainsize)),         nrow = trainsize,       ncol = nsnps)
Ztest     =  matrix(rep(0, (nsnps*testsize)),          nrow = testsize,        ncol = nsnps)
ZZtrain   =  matrix(rep(0, (features*trainsize)),         nrow = trainsize,    ncol = features)
ZZtest    =  matrix(rep(0, (features*testsize)),          nrow = testsize,     ncol = features)





Ptrain<-c();
for(j in 1:nsnps){
if(k == 1) {Ztrain[,j]<-Z[(1+sfold):samples, j]; Ztest[,j]<-Z[1:sfold, j]; Ptrain<-P[(1+sfold):samples];}
else{
if(k == 10){Ztrain[,j]<-Z[1:strain,j]; Ztest[,j]<-Z[(1+strain):samples, j]; Ptrain<-P[1:strain];}
else{
Ztrain[,j]<-c(Z[1:(sfold*k - sfold), j], Z[(sfold*k+1):samples,j]); Ptrain<-c(P[1:(sfold*k-sfold)],P[(sfold*k+1):samples]);
 Ztest[,j]<-Z[(sfold*k - sfold + 1):(sfold*k),j];
}}
}


corlist<-c();for(snp in 1:nsnps){corlist[snp] = cor(Ztrain[,snp], Ptrain)}
top<-c(); for(f in 1:features){top[f] = sort(corlist)[nsnps - f + 1]}


flag<-rep(0, nsnps); index<-c();
for(f in 1:features){
for(snp in 1:nsnps){if(corlist[snp] == top[f] & flag[snp] == 0){index[f] = snp; flag[snp] = 1; break;}}
ZZtrain[,f] = Ztrain[,index[f]]
ZZtest[,f] = Ztest[,index[f]]
}
ZZtrain = provideDimnames(ZZtrain, sep = "", base = list(LETTERS), unique = TRUE)
ZZtest  = provideDimnames(ZZtest, sep = "",  base = list(LETTERS), unique = TRUE)




if(trait == 2){
y          = as.factor(Ptrain)
M1         = svm(ZZtrain, y, kernel = "linear", type = "C-classification"); Psvm = as.numeric(predict(M1, ZZtrain)) - 1;
lgbtrain   = lgb.Dataset(data=as.matrix(ZZtrain), label=y, colnames = colnames(ZZtrain))
M2         = lgb.train(params=list(objective="binary",learning_rate=0.01,boosting="gbdt",num_threads=6,seed=1238845,num_leaves=5,verbosity=-1),data=lgbtrain,nrounds=10000,colnames=colnames(ZZtrain))
Prf        = as.numeric(predict(M2, data = as.matrix(ZZtrain)))
M3         = glm(Ptrain ~., data.frame(ZZtrain), family=binomial, control=list(maxit=100,trace=F));
Pmr        = as.numeric(predict(M3, data.frame(ZZtrain), type=c("response")));
A          = data.frame(matrix(c(Psvm, Prf, Pmr, as.numeric(ZZtrain)), ncol = 3 + features, byrow = FALSE));
A          = provideDimnames(A,sep = "", base = list(LETTERS), unique = TRUE);
lgbblend   = lgb.Dataset(data=as.matrix(A), label=y,colnames=colnames(A))
Mblend     = lgb.train(params=list(objective="binary",learning_rate=0.01,boosting="gbdt",num_threads=6,seed=1238845,num_leaves=5,verbosity=-1),data=lgbblend,nrounds=10000,colnames=colnames(A))

Psvmtest   = as.numeric(predict(M1, ZZtest)) - 1;
Prftest    = as.numeric(predict(M2, data = as.matrix(ZZtest)))
Pmrtest    = as.numeric(predict(M3, data.frame(ZZtest), type=c("response")));
Atest      = data.frame(matrix(c(Psvmtest, Prftest, Pmrtest, as.numeric(ZZtest)), ncol = 3 + features, byrow = FALSE));
Atest      = provideDimnames(Atest, sep = "", base = list(LETTERS), unique = TRUE);
Ptemp      = as.numeric(predict(Mblend, data = as.matrix(Atest)))
}


if(trait == 1){
y          = as.numeric(Ptrain);
M1         = svm(ZZtrain, Ptrain, kernel = "linear", type = "eps-regression");    Psvm = as.numeric(predict(M1, ZZtrain));
lgbtrain   = lgb.Dataset(data=as.matrix(ZZtrain), label=y,colnames=colnames(ZZtrain))
M2         = lgb.train(params=list(objective="regression",learning_rate=0.01,boosting="gbdt",num_threads=6,seed=1238845,num_leaves=5,verbosity=-1),data=lgbtrain,nrounds=10000,colnames=colnames(ZZtrain))
Prf        = as.numeric(predict(M2, data = as.matrix(ZZtrain)))
M3         = lm(Ptrain ~., data.frame(ZZtrain)); Pmr  = as.numeric(predict(M3, data.frame(ZZtrain)))
A          = data.frame(matrix(c(Psvm, Prf, Pmr, as.numeric(ZZtrain)), ncol = 3 + features, byrow = FALSE))
A          = provideDimnames(A,sep = "", base = list(LETTERS), unique = TRUE)
lgbblend   = lgb.Dataset(data=as.matrix(A), label=y,colnames=colnames(A))
Mblend     = lgb.train(params=list(objective="regression",learning_rate=0.01,boosting="gbdt",num_threads=6,seed=1238845,num_leaves=5,verbosity=-1),data=lgbblend,nrounds=10000,colnames=colnames(A))

Psvmtest   = as.numeric(predict(M1, ZZtest));
Prftest    = as.numeric(predict(M2, data = as.matrix(ZZtest)))
Pmrtest    = as.numeric(predict(M3, data.frame(ZZtest)))
Atest      = data.frame(matrix(c(Psvmtest, Prftest, Pmrtest, as.numeric(ZZtest)), ncol = 3 + features, byrow = FALSE))
Atest      = provideDimnames(Atest, sep = "", base = list(LETTERS), unique = TRUE)
Ptemp      = as.numeric(predict(Mblend, data = as.matrix(Atest)))
}




filename   = paste("mmassoc-cvfold", k, ".txt", sep=""); 
write(Ptemp, ncolumns=length(Ptemp), file=filename);
}




Pfinal<-c(); 
for(k in 1:10){
filename   = paste("mmassoc-cvfold", k, ".txt", sep=""); 
Ptemp      = scan(filename);
Pfinal     = c(Pfinal, Ptemp);
}




if(covs == 0){
COV  = matrix(runif(samples), ncol = 1, nrow = samples); 
covs = 1;
}




if(trait == 1){
pcor1 = cor(Pfinal, P, method =  "pearson");
pcor2 = cor(Pfinal, P, method = "spearman");
pcor3 = cor(Pfinal, P, method =  "kendall");
res   = lm(P ~ Pfinal); beta = as.numeric(res$coefficients[2]);

if(pcor1 >= 0 & pcor2 >=0 & pcor3 >= 0 & beta >= 0){
ZP = data.frame(matrix(c(Pfinal, as.numeric(Z), as.numeric(COV)), ncol = (1 + nsnps + covs), byrow = FALSE));
fitall<-lm(P ~., ZP              ); 
fitcov<-lm(P ~., data.frame(COV) ); 
plr = (anova(fitcov, fitall)[6]$`Pr(>F)`)[2];
remove(ZP); remove(fitall); remove(fitcov); 
}
 
else{
ZP = data.frame(matrix(c(as.numeric(Z), as.numeric(COV)), ncol = (nsnps + covs), byrow = FALSE));
fitall<-lm(P ~., ZP              ); 
fitcov<-lm(P ~., data.frame(COV) ); 
plr = (anova(fitcov, fitall)[6]$`Pr(>F)`)[2];
remove(ZP); remove(fitall); remove(fitcov); 
}

remove(pcor1); remove(pcor2); remove(pcor3); remove(res); remove(beta);
}




if(trait == 2){
pcor1  = cor(Pfinal, P, method = "pearson");
pcor2  = cor(Pfinal, P, method = "spearman");
pcor3  = cor(Pfinal, P, method = "kendall");
logreg = glm(P ~ Pfinal, family = binomial); beta = as.numeric(logreg[1]$coefficients[2]);

if(pcor1 >= 0 & pcor2 >=0 & pcor3 >= 0 & beta >= 0){
ZP = data.frame(matrix(c(Pfinal, as.numeric(Z), as.numeric(COV)), ncol = (1 + nsnps + covs), byrow = FALSE));
fitall<-glm(P ~., ZP,              family = binomial); resall = summary(fitall);
fitcov<-glm(P ~., data.frame(COV), family = binomial); rescov = summary(fitcov);
plr = pchisq((as.numeric(rescov[4]) - as.numeric(resall[4])), (1 + nsnps), lower.tail = FALSE);
remove(ZP); remove(fitall); remove(fitcov); remove(resall); remove(rescov);
}

else{
ZP = data.frame(matrix(c(as.numeric(Z), as.numeric(COV)), ncol = (nsnps + covs), byrow = FALSE));
fitall<-glm(P ~., ZP,              family = binomial); resall = summary(fitall);
fitcov<-glm(P ~., data.frame(COV), family = binomial); rescov = summary(fitcov);
plr = pchisq((as.numeric(rescov[4]) - as.numeric(resall[4])), (nsnps), lower.tail = FALSE);
remove(ZP); remove(fitall); remove(fitcov); remove(resall); remove(rescov);
}

remove(pcor1); remove(pcor2); remove(pcor3); remove(logreg); remove(beta);
}






result = paste("P value = ", plr, sep="");
write(result, file = "pvalues.txt", append = T);


















