library(foreach);
library(doParallel);
registerDoParallel(cores=detectCores());
library(ranger);
library(e1071); 

args      <- commandArgs(trailingOnly = TRUE);
nsnps     =   as.integer(args[2]);
features  =   as.integer(args[4]);
samples   =   as.integer(args[6]);
trait     =   as.integer(args[8]);
covs      =   as.integer(args[10]);
in.snp    =   as.character(args[12]);
in.pheno  =   as.character(args[14]);
in.covs   =   as.character(args[16]);

sfold     =   as.integer(0.1*samples);
strain    =   9*sfold;
Z         =   matrix(scan(in.snp),  nrow = samples,  ncol = nsnps, byrow=T);
P         =   scan(in.pheno);
if(covs > 0){
COV       =   matrix(scan(in.covs), nrow = samples,  ncol = covs,  byrow=T);
}
reorder   =   sample(c(1:samples));





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





Ztrain    =  matrix(rep(0, (nsnps*trainsize)),            nrow = trainsize,       ncol = nsnps);
Ztest     =  matrix(rep(0, (nsnps*testsize)),             nrow = testsize,        ncol = nsnps);
ZZtrain   =  matrix(rep(0, (features*trainsize)),         nrow = trainsize,    ncol = features);
ZZtest    =  matrix(rep(0, (features*testsize)),          nrow = testsize,     ncol = features);





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


corlist<-c();for(snp in 1:nsnps){corlist[snp] = cor(Ztrain[,snp], Ptrain);}
top<-c(); for(f in 1:features){top[f] = sort(corlist)[nsnps-f+1];}


flag<-rep(0, nsnps); index<-c();
for(f in 1:features){
for(snp in 1:nsnps){if(corlist[snp] == top[f] & flag[snp] == 0){index[f] = snp; flag[snp] = 1; break;}}
ZZtrain[,f] = Ztrain[,index[f]]; ZZtest[,f] = Ztest[,index[f]];
}

ZZtrain = provideDimnames(ZZtrain, sep = "", base = list(LETTERS), unique = TRUE);
ZZtest  = provideDimnames(ZZtest, sep = "",  base = list(LETTERS), unique = TRUE);


if(trait == 2){
y          = as.factor(Ptrain);
M1         = svm(ZZtrain, as.factor(Ptrain), kernel = "linear", type = "C-classification"); Psvm = as.numeric(predict(M1, ZZtrain)) - 1;
M2         = ranger(y ~., data = data.frame(ZZtrain), mtry = 1, num.trees = 1001, write.forest=T);      
Prf        = as.numeric(predict(M2, data.frame(ZZtrain), num.trees = M2$num.trees, type = "response")$predictions) - 1; 
M3         = glm(y ~., data.frame(ZZtrain), family=binomial, control=list(maxit=100,trace=F)); 
Pmr  = as.numeric(predict(M3, data.frame(ZZtrain), type=c("response")));
A          = data.frame(matrix(c(Psvm, Prf, Pmr, as.numeric(ZZtrain)), ncol = 3 + features, byrow = FALSE));
A          = provideDimnames(A,sep = "", base = list(LETTERS), unique = TRUE);
Mblend     = ranger(y ~., data = A, mtry = 1, num.trees = 1001, write.forest=T);      
Psvmtest   = as.numeric(predict(M1, ZZtest)) - 1;
Prftest    = as.numeric(predict(M2, data.frame(ZZtest), num.trees = M2$num.trees, type = "response")$predictions) - 1;
Pmrtest    = as.numeric(predict(M3, data.frame(ZZtest), type=c("response")));
Atest      = data.frame(matrix(c(Psvmtest, Prftest, Pmrtest, as.numeric(ZZtest)), ncol = 3 + features, byrow = FALSE));
Atest      = provideDimnames(Atest, sep = "", base = list(LETTERS), unique = TRUE);
Ptemp      = as.numeric(predict(Mblend, Atest, num.trees = Mblend$num.trees, type = "response")$predictions) - 1; 
}

if(trait == 1){
y          = as.numeric(Ptrain);
M1         = svm(ZZtrain, Ptrain, kernel = "linear", type = "eps-regression");    Psvm = as.numeric(predict(M1, ZZtrain));
M2         = ranger(y ~., data = data.frame(ZZtrain), mtry = 1, num.trees = 1001, write.forest=T);      
Prf        = as.numeric(predict(M2, data.frame(ZZtrain), num.trees = M2$num.trees, type = "response")$predictions); 
M3         = lm(y ~., data.frame(ZZtrain)); Pmr  = as.numeric(predict(M3, data.frame(ZZtrain)));
A          = data.frame(matrix(c(Psvm, Prf, Pmr, as.numeric(ZZtrain)), ncol = 3 + features, byrow = FALSE));
A          = provideDimnames(A,sep = "", base = list(LETTERS), unique = TRUE);
Mblend     = ranger(y ~., data = A, mtry = 1, num.trees = 1001, write.forest=T);      
Psvmtest   = as.numeric(predict(M1, ZZtest));
Prftest    = as.numeric(predict(M2, ZZtest, num.trees = M2$num.trees, type = "response")$predictions);
Pmrtest    = as.numeric(predict(M3, data.frame(ZZtest)));
Atest      = data.frame(matrix(c(Psvmtest, Prftest, Pmrtest, as.numeric(ZZtest)), ncol = 3 + features, byrow = FALSE));
Atest      = provideDimnames(Atest, sep = "", base = list(LETTERS), unique = TRUE);
Ptemp      = as.numeric(predict(Mblend, Atest, num.trees = Mblend$num.trees, type = "response")$predictions);
}

filename   = paste("mmassoc-cvfold-ensemble", k, ".txt", sep="");   write(Ptemp,   ncolumns=length(Ptemp),   file=filename);
filename   = paste("mmassoc-cvfold-linear",   k, ".txt", sep="");   write(Pmrtest, ncolumns=length(Pmrtest), file=filename);
}




Pfinal<-c(); Plinear<-c();
for(k in 1:10){
filename   = paste("mmassoc-cvfold-ensemble", k, ".txt", sep=""); 
Ptemp      = scan(filename);
Pfinal     = c(Pfinal, Ptemp);
filename   = paste("mmassoc-cvfold-linear", k, ".txt", sep=""); 
Ptemp      = scan(filename);
Plinear    = c(Plinear, Ptemp);
}




if(covs == 0){
COV  = matrix(runif(samples), ncol = 1, nrow = samples); 
covs = 1;
}

pcor1 = cor(P,  Pfinal);
pcor2 = cor(P, Plinear);

if(trait == 1){
ZINT = data.frame(matrix(c(Pfinal, Plinear, as.numeric(Z), as.numeric(COV)), ncol = (2 + nsnps + covs),     byrow = FALSE));
ZLIN = data.frame(matrix(c(Plinear, as.numeric(Z), as.numeric(COV)),         ncol = (1 + nsnps + covs),     byrow = FALSE));
fitint<-lm(P ~., ZINT); 
fitlin<-lm(P ~., ZLIN); 
if(pcor1 >= pcor2){
plr = (anova(fitlin, fitint)[6]$`Pr(>F)`)[2];
}
else{plr = NA;}
remove(ZINT); remove(ZLIN); remove(fitint); remove(fitall); 
}




if(trait == 2){
ZINT = data.frame(matrix(c(Pfinal,  Plinear, as.numeric(Z), as.numeric(COV)), ncol = (2 + nsnps + covs),     byrow = FALSE));
ZLIN = data.frame(matrix(c(Plinear, as.numeric(Z), as.numeric(COV)),          ncol = (1 + nsnps + covs),     byrow = FALSE));
fitint<-glm(P ~., ZINT,   family = binomial); resint = summary(fitint);
fitlin<-glm(P ~., ZLIN,   family = binomial); reslin = summary(fitlin);
if(pcor1 >= pcor2){
plr = pchisq((as.numeric(reslin[4]) - as.numeric(resint[4])), 1, lower.tail = FALSE);
}
else{plr = NA;}
remove(ZINT); remove(ZLIN); remove(fitint); remove(fitlin); remove(resint); remove(reslin);
}






result = paste("Interactions P value = ", plr, sep=""); print(result);


















