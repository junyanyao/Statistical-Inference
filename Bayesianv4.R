require(foreign)

#dat<- read.csv("C:/Users/jyao/Downloads/magnesium.csv", header = TRUE)
dat<- read.csv("~/Downloads/magnesium.csv", header = TRUE)
require(metafor)
require(rmeta)

######################peto analysis###################################
#peto method for all 15 trials
peto15<-rma.peto(ai=dead1, n1i=tot1,ci=dead0,n2i =tot0, data=dat)
peto15

#peto method for 8 trials
cat<- sample(0,15, replace = TRUE)
dat2<- cbind(dat, cat)
dat2$cat[1:7]<-1
dat2$cat[13]<-1


peto8<-rma.peto(ai=dead1, n1i=tot1,ci=dead0,n2i=tot0, data=dat2, subset = (cat==1))
peto8

#peto method for 14 trails
cat2<- sample(0,15, replace = TRUE)
dat3<- cbind(dat2, cat2)
dat3$cat2[1:14]<-1

peto14<-rma.peto(ai=dead1, n1i=tot1,ci=dead0,n2i=tot0, data=dat3, subset = (cat2==1))
peto14


#alternative way
a<- escalc(measure = "PETO",ai=dead1, n1i=tot1,ci=dead0,n2i=tot0, data=dat, add = 0 )
b<- rma(yi,vi, data=a, method = "FE")

#Forest plot
forest(a$yi,a$vi, slab = dat$trialnam, transf = exp)
forest(a$yi,a$vi, slab = dat$trialnam)

#######################DerSimonian Laird analysis######################
#all 15 trials
DL15<- meta.DSL(ntrt = tot1, nctrl = tot0, ptrt = dead1, pctrl = dead0, data = dat)
DL15

#8 trials
DL8<- meta.DSL(ntrt = tot1, nctrl = tot0, ptrt = dead1, pctrl = dead0, data = dat2, subset=(cat==1))
dat99<- dat2[dat2$cat==1,]
DL8<- meta.DSL(ntrt = tot1, nctrl = tot0, ptrt = dead1, pctrl = dead0, data = dat99)
DL8

#14 trials
dat88<- dat3[dat3$cat2==1,]
DL14<-meta.DSL(ntrt = tot1, nctrl = tot0, ptrt = dead1, pctrl = dead0, data = dat88)


######Bayesian Parts #################

k=15
nc<-dat$tot0
nm<-dat$tot1
rc=dat$dead0
rm=dat$dead1

set.seed(100)
model = "data {
int <lower = 0> k; 
int <lower = 0> nc [k];
int <lower = 0> nm [k];
int <lower = 0> rc [k];
int <lower = 0> rm [k];
}

parameters{
real <lower = 0, upper = 1> pc[k];
vector[k] delta;
real mu;
real <lower = 0> sigma;
real deltanew;
}

transformed parameters {
real <lower = 0, upper = 1> pm[k];
for (i in 1:k) {
pm[i] = exp(log(pc[i]/(1-pc[i])) + delta[i])/(1+exp(log(pc[i]/(1-pc[i])) + delta[i]));
}
}

model {
# model
#for (i in 1:k) {
rc ~ binomial(nc, pc);
rm ~ binomial(nm, pm);
delta ~ normal(mu, sigma);
#}
deltanew ~ normal(mu, sigma);

# priors
#for (i in 1:k) {
pc ~ uniform(0,1);
#}

mu ~ normal(0, 100);
sigma ~ uniform(0, 100);

}"

fit1 = stan(model_code = model, data = c("k", "nc", "nm", "rc", "rm"), pars = "deltanew", iter = 500000, chains = 3, warmup = 500)
print(fit1)
traceplot(fit1, pars= "deltanew")

#odd ratio =0.3945

########################SKEPTICAL#################################

model = "data {
int <lower = 0> k; 
int <lower = 0> nc [k];
int <lower = 0> nm [k];
int <lower = 0> rc [k];
int <lower = 0> rm [k];
}

parameters{
real <lower = 0, upper = 1> pc[k];
vector[k] delta;
real mu;
real <lower = 0> sigma;
real deltanew;
}

transformed parameters {
real <lower = 0, upper = 1> pm[k];
for (i in 1:k) {
pm[i] = exp(log(pc[i]/(1-pc[i])) + delta[i])/(1+exp(log(pc[i]/(1-pc[i])) + delta[i]));
}
}

model {
# model
#for (i in 1:k) {
rc ~ binomial(nc, pc);
rm ~ binomial(nm, pm);
delta ~ normal(mu, sigma);
#}
deltanew ~ normal(mu, sigma);

# priors
#for (i in 1:k) {
pc ~ uniform(0,1);
#}

mu ~ normal(0, 0.1749);
sigma ~ uniform(0, 100);

}"

fit2 = stan(model_code = model, data = c("k", "nc", "nm", "rc", "rm"), pars = "deltanew", iter = 500000, chains = 3, warmup = 500)
print(fit2)
#odd ratio =0.74
traceplot(fit2, pars= "deltanew")


efit1 <- extract (fit1, permuted=TRUE)
hist(exp(efit1$deltanew), xlim =c(0,2), breaks = 10000000, main="Reference")


efit2 <- extract (fit2, permuted=TRUE)
hist(exp(efit2$deltanew), xlim = c(0,2), breaks = 10000000, main="Skeptical")

