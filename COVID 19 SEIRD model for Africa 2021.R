
###---------setting the working directory-----------------

setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")


## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")

###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")

###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2020-03-01")))/(365.25)+2020
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################


head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New


Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor) #infection fatality rate
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- 1e-5 # the initial immunization level
r_in <- 0.145 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate

## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")

dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")

# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun)*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop, birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.65 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.3 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new

Date=seq(as.Date("2020/3/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
#write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))
DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))
lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)
CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)

C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)

#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"CurrentAdupdateedMarch.xlsx", sheetName =COUNTRY, append=TRUE)

}


##arranging the data in Ms Excel
sheet = excel_sheets("CurrentAdupdateedMarch.xlsx")
Combined_country_dataMarch = lapply(setNames(sheet, sheet), function(x) read_excel("CurrentAdupdateedMarch.xlsx", sheet=x))
Combined_country_dataMarch = bind_rows(Combined_country_dataMarch, .id="Sheet")
Combined_country_dataMarch =as.data.frame(Combined_country_dataMarch)
Combined_country_dataMarch$Date <- as.Date(Combined_country_dataMarch$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_dataMarch, "Combined_country_dataMarch.xlsx")
head(Combined_country_dataMarch)



##############################################################################
###########################################################################

#################################################################
#
######projection for 2022 starting from the current immunized populations per country
#  ...at the current rate of immunization and rates of immization and variants
####################################################################

setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")



## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


###---------setting the working directory-----------------

##setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter")

######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")

###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")

###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################

head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New

Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate, omicron odds ratio to delta in 0.3, implying a reducition of 70%
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate

## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")
dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")

# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################

times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new


Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"Current22.xlsx", sheetName=COUNTRY, append=TRUE)

}


sheet = excel_sheets("Current22.xlsx")
Combined_country_data22 = lapply(setNames(sheet, sheet), function(x) read_excel("Current22.xlsx", sheet=x))
Combined_country_data22 = bind_rows(Combined_country_data22, .id="Sheet")
Combined_country_data22 =as.data.frame(Combined_country_data22)
Combined_country_data22$Date <- as.Date(Combined_country_data22$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_data22, "Combined_country_data22.xlsx")
head(Combined_country_data22)



##############################################################################
#### case scenarios

##############################################################################
###
##  Re - infection

#### A new variant that is more reinfectious than the delta by 80%, 120% and 200%

####### change the ephsilonh at those levels in the DynCovInitGAReinf.txt initializing data set with a dominance of 80%

##########################################################################


##  80% more

#######

setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")


## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)



######---------- load the packages needed, and set the random seed, to allow reproducibility

library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGAReinf.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")



###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################

head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New


Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate



## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")



dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonh80, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new


Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported

lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


#Case=melt(CasesH, id.vars=c("Date"))
#head(Case)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

#library(data.table)
#library("Hmisc")
#library(xlsx)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"CurrentEpsi80.xlsx", sheetName =COUNTRY, append=TRUE)
}


sheet = excel_sheets("CurrentEpsi80.xlsx")
Combined_country_dataEpsi80 = lapply(setNames(sheet, sheet), function(x) read_excel("CurrentEpsi80.xlsx", sheet=x))
Combined_country_dataEpsi80 = bind_rows(Combined_country_dataEpsi80, .id="Sheet")
Combined_country_dataEpsi80 =as.data.frame(Combined_country_dataEpsi80)
Combined_country_dataEpsi80$Date <- as.Date(Combined_country_dataEpsi80$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_dataEpsi80, "Combined_country_dataEpsi80.xlsx")
head(Combined_country_dataEpsi80)


##########################################

## 120% more

##########################################

setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()


Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")

## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)

######---------- load the packages needed, and set the random seed, to allow reproducibility

library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGAReinf.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")


as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")


###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################


head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New


Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate

## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")

dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonh120, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new


Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported

lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))

Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


#Case=melt(CasesH, id.vars=c("Date"))
#head(Case)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"CurrentEpsi120.xlsx", sheetName =COUNTRY, append=TRUE)
}


sheet = excel_sheets("CurrentEpsi120.xlsx")
Combined_country_dataEpsi120 = lapply(setNames(sheet, sheet), function(x) read_excel("CurrentEpsi120.xlsx", sheet=x))
Combined_country_dataEpsi120 = bind_rows(Combined_country_dataEpsi120, .id="Sheet")
Combined_country_dataEpsi120 =as.data.frame(Combined_country_dataEpsi120)
Combined_country_dataEpsi120$Date <- as.Date(Combined_country_dataEpsi120$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_dataEpsi120, "Combined_country_dataEpsi120.xlsx")
head(Combined_country_dataEpsi120)


##########################################

####   200% more

##########################################


setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")

## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


######---------- load the packages needed, and set the random seed, to allow reproducibility
library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGAReinf.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")




###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################

head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New


Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate



## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")

dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonh200, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new


Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported

lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))





Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


#Case=melt(CasesH, id.vars=c("Date"))
#head(Case)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


sum(Cases) 
sum(H)

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"CurrentEpsi200.xlsx", sheetName =COUNTRY, append=TRUE)
}


sheet = excel_sheets("CurrentEpsi200.xlsx")
Combined_country_dataEpsi200 = lapply(setNames(sheet, sheet), function(x) read_excel("CurrentEpsi200.xlsx", sheet=x))
Combined_country_dataEpsi200 = bind_rows(Combined_country_dataEpsi200, .id="Sheet")
Combined_country_dataEpsi200 =as.data.frame(Combined_country_dataEpsi200)
Combined_country_dataEpsi200$Date <- as.Date(Combined_country_dataEpsi200$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_dataEpsi200, "Combined_country_dataEpsi200.xlsx")
head(Combined_country_d ataEpsi200)



##############################################################################
#### case scenarios

##############################################################################
###
##  Variant severity

#### A new variant that is more severe than the delta by  80%, 120% and 200% 

####### adjust the rates[] 8 and 9 accordingly uising in the DynCovInitGA.txt initializing data set with a dominance of 70%

##########################################################################





######
#### 80% more
######



setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")



## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


###---------setting the working directory-----------------

##setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter")

######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")




###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################

head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New

Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate, omicron odds ratio to delta in 0.3, implying a reducition of 70%
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate



## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")



dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = 1.8*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = 1.8*(1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new


Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

#library(data.table)
#library("Hmisc")
#library(xlsx)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"Currentsev80.xlsx", sheetName=COUNTRY, append=TRUE)
}

##############################################################################
###
##  Variant severity

#### A new variant that is more severe than the delta by  80%, 120% and 200% 

####### adjust the rates[] 8 and 9 accordingly uising in the DynCovInitGA.txt initializing data set with a dominance of 70%

##########################################################################


######
#### 120% more
######




setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")



## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


###---------setting the working directory-----------------

##setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter")

######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

#write.xlsx(covinit, "DynCovInitGA.xlsx")



#filter(covcase, Country==COUNTRY)

#write.csv(covinit, "covinit.csv")
#write.csv(covcase, "covcase.csv")





as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")




###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################


head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New





Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate, omicron odds ratio to delta in 0.3, implying a reducition of 70%
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate



## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")



dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = 2.2*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = 2.2*(1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new


Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))





Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


#Case=melt(CasesH, id.vars=c("Date"))
#head(Case)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

#library(data.table)
#library("Hmisc")
#library(xlsx)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"Currentsev120.xlsx", sheetName=COUNTRY, append=TRUE)
}


##############################################################################
###
##  Variant severity

#### A new variant that is more severe than the delta by  80%, 120% and 200% 

####### adjust the rates[] 8 and 9 accordingly uising in the DynCovInitGA.txt initializing data set with a dominance of 70%

##########################################################################




##### 200% more



setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")



## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


###---------setting the working directory-----------------

##setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter")

######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")




###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################


head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New


Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate, omicron odds ratio to delta in 0.3, implying a reducition of 70%
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans+future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate



## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")



dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = 3*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = 3*(1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new

Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


#Case=melt(CasesH, id.vars=c("Date"))
#head(Case)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


sum(Cases) 
sum(H)

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

#library(data.table)
#library("Hmisc")
#library(xlsx)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"Currentsev200.xlsx", sheetName=COUNTRY, append=TRUE)
}


sheet = excel_sheets("Currentsev80.xlsx")
Combined_country_datasev80 = lapply(setNames(sheet, sheet), function(x) read_excel("Currentsev80.xlsx", sheet=x))
Combined_country_datasev80 = bind_rows(Combined_country_datasev80, .id="Sheet")
Combined_country_datasev80 =as.data.frame(Combined_country_datasev80)
Combined_country_datasev80$Date <- as.Date(Combined_country_datasev80$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_datasev80, "Combined_country_datasev80.xlsx")
head(Combined_country_datasev80)

sheet = excel_sheets("Currentsev120.xlsx")
Combined_country_datasev120 = lapply(setNames(sheet, sheet), function(x) read_excel("Currentsev120.xlsx", sheet=x))
Combined_country_datasev120 = bind_rows(Combined_country_datasev120, .id="Sheet")
Combined_country_datasev120 =as.data.frame(Combined_country_datasev120)
Combined_country_datasev120$Date <- as.Date(Combined_country_datasev120$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_datasev120, "Combined_country_datasev120.xlsx")
head(Combined_country_datasev120)

sheet = excel_sheets("Currentsev200.xlsx")
Combined_country_datasev200 = lapply(setNames(sheet, sheet), function(x) read_excel("Currentsev200.xlsx", sheet=x))
Combined_country_datasev200 = bind_rows(Combined_country_datasev200, .id="Sheet")
Combined_country_datasev200 =as.data.frame(Combined_country_datasev200)
Combined_country_datasev200$Date <- as.Date(Combined_country_datasev200$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_datasev200, "Combined_country_datasev200.xlsx")
head(Combined_country_datasev200)








##############################################################################
#### case scenarios

##############################################################################
###
##  Transmissibility

#### A new variant that is more transmissible than the delta by  80%, 120% and 200% 

####### increase the future1_trans by the percentage increases in the DynCovInitGATrans.txt initializing data set with a dominance of 80%

##########################################################################

#### 80%


setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")



## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


###---------setting the working directory-----------------

##setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter")

######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")

###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################

head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New


Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate, omicron odds ratio to delta in 0.3, implying a reducition of 70%
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans*1.8 +future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate

## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")

dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")
# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new


Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))

Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


sum(Cases) 
sum(H)

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

#library(data.table)
#library("Hmisc")
#library(xlsx)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"CurrentTrans80.xlsx", sheetName=COUNTRY, append=TRUE)

}







##################################################################

###### 120%

##################################################################

setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")



## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


###---------setting the working directory-----------------

##setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter")

######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")


as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")

###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################
head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New

Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate, omicron odds ratio to delta in 0.3, implying a reducition of 70%
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans*2.2 +future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate



## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")



dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new

Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


sum(Cases) 
sum(H)

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"CurrentTrans120.xlsx", sheetName=COUNTRY, append=TRUE)

}



################################################


##### 200% More


#################################################


setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter/Country_Models_Gamma_Noise/Improved two in one model")
rm(list = ls()) # removing the objects in the memory
ls()

Co=c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde",
"Cameroon","Central African Republic","Chad","Comoros","Congo","Cote d Ivoire",
"Democratic Republic of the Congo","Equatorial Guinea","Eritrea","Eswatini",
"Ethiopia","Gabon","Gambia","Ghana","Guinea","GuineaBissau","Kenya",
"Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria",
"Rwanda","Sao Tome and Principe","Senegal",
"Seychelles","Sierra Leone","South Africa","South Sudan","Togo","Uganda",
"United Republic of Tanzania","Zambia","Zimbabwe")



## ----opts,include=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


###---------setting the working directory-----------------

##setwd("F:/Epi Modelling COVID-19/Data/Analytics/Filter")

######---------- load the packages needed, and set the random seed, to allow reproducibility


library(devtools)
#PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# installing the required libraries
library(pomp)
library(ggplot2)
library(zoo)
library(fitdistrplus)
library(rJava)
library(xlsx)
library(openxlsx)
library(tidyverse)
library(viridis)
library(dplyr)
library(growthcurver)
library(tibbletime)
library(tidyquant)
library(readxl)
library(tidyverse)
library(scales)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)


## ----select-Country,echo=FALSE---------------------------------------------------
for(COUNTRY in Co)

{
#####-------select a Country------------------

## ----select-Country,echo=FALSE---------------------------------------------------
print(COUNTRY)

####-----loading data and covarites

###--importing------

covcase=read.table("CovCasesWeekly.txt", header=TRUE,sep="\t")
head(covcase)
attach(covcase)
covinit=read.table("DynCovInitGA.txt", header=TRUE,sep="\t")
head(covinit)
covcase$date <- as.Date(date, "%d/%m/%Y")

as.tibble(covcase, key=Country) -> covtib 
#covcaset <- as_tsibble(covtib, index = date)
####----saving .rda objects

# Saving a R object in RData format
getwd()
save(covtib, file = "covcases.rda")
load("covcases.rda")


###----the Country to be selected at this point
##e.g Country <- "Kenya" or any other Country, If not selected, then it should
# just read Country and wait for the selection

## ----load-data----------------------------------------------------------------
load("covcases.rda")
covtib %>%
  filter(Country==COUNTRY) -> covdat

#plot(covdat$new_cases)
#plot the data and covariates.
covdat %>% ggplot(aes(x=date,y=new_cases))+geom_line()+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
 xlab("Time in Months") + ylab("Confirmed Cases")



###----------filtering the countries
covcase %>%
 mutate(year=as.integer(format(date,"%Y"))) %>%
  filter(Country==COUNTRY ) %>%
  mutate(
    time=(julian(date,origin=as.Date("2022-01-01")))/(365.25)+2022
  ) -> cdat

######-----filtering initial parameters
covinit %>%
  filter(Country==COUNTRY)-> covinit_par

attach(covinit_par)



##########################################
#
####Country: Country XXX
#
##########################################

###Fiting the fourier series to data and extracting the parameters

##########################################################
#                                                        #
#          Model with the two sub-models                 #
#                                                        #
##########################################################


head(cdat)
attach(cdat)

Nl=dim(cdat)
NL=Nl[1]

NL=Nl[1]
casN_Imm=ceiling(cdat[[3]]*(1-prop_pos_imm))#proportion of cases not - immunized
fit11 <- fitdist(casN_Imm, "nbinom")
summary(fit11)
thetah1=fit11$estimate[[1]]
rhoH1=fit11$estimate[[2]]
thetah1
rhoH1


casImm=round(cdat[[3]]*(prop_pos_imm)) ##proportion of the cases from the immunized
fit12 <- fitdist(casImm, "nbinom")
summary(fit12)
thetah2=fit12$estimate[[1]]
rhoH2=fit12$estimate[[2]]
thetah2
rhoH2


#death<-read.table("deaths.txt", header=TRUE,sep="\t")
fit2 <- fitdist(cdat[[4]], "nbinom")
lam=fit2$estimate[[1]]
rhoD_New=fit2$estimate[[2]]
lam
rhoD_New


Pop
#### Initializing the parameters
ref_var_prop = 1- (alpha_prop+beta_prop+delta_prop+future1_prop+future2_prop) #proportion of  the reference variant
CFR = (Dcases/TestsP)*(lenth_of_stay/365.25)*(vun_pop_infl_factor)*0.3 #infection fatality rate, omicron odds ratio to delta in 0.3, implying a reducition of 70%
CFR1 = CFR #rate per two weeks (time to event)
latent=(latent_alpha*alpha_prop + latent_beta*beta_prop+
latent_delta*delta_prop+latent_future1*future1_prop+latent_referece*ref_var_prop)
latent
lambdah = 1/latent #incubation period (1/latent period)
gammah1 = 1 / infectious_period  #recovery rate
gammah2 = 1 / (infectious_period-2.8)  #recovery rate for the immunized


#Deriving immunization rates
times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)
n0_in <- prop_pop_imm_21 # the initial immunization level
r_in <- 0.09 # the initial growth rate
N <- length(tt) # the number of "measurements" collected during the growth
data_t <- 0:N  # the times the measurements were made
imm <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
 #Now summarize the "experimental" growth data that we just generated
gc <- SummarizeGrowth(data_t, imm)
plot(gc)

#datesd=0:(dim(cdat)[1]-1)
#gc_fit <- SummarizeGrowth(datesd, cdat$new_cases)

## Construct some birthrate data.

data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    birthrate=Pop*bspline.basis(time,nbasis=5)%*%(rep(Anbirthrate[1],5))
  ) -> birthdat

### Population proportions between the general population and the immunized
data.frame(time=seq(-1,11.4,by=1/12)) %>%
  mutate(
    immunized=(NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t))*Pop 
  ) -> immun


prop_imm = imm
Pop1 = Pop*(1-prop_imm)
plot(Pop1)
#imm=prop_imm
prop_imm = imm
plot(prop_imm)


#calculating variant severities
severity_var = (alpha_prop*alpha_trans+beta_prop*beta_trans+delta_prop*delta_trans
   +ref_var_prop*reference_trans+future1_prop*future1_trans*3 +future2_prop*future2_trans)


#Obtaining the seasonality in the data
y=log(cdat$new_cases+1)
y


DataC <- data.frame(x = seq_along(y), y = y)
N <- length(DataC$x)
L <- DataC$x[N]-DataC$x[1]
trend <- lm(y ~ x, data = DataC)
DataC$detrended <- DataC$y-trend$fitted.values
fft_in <- DataC$detrended
fft_out <- fft(fft_in)
barplot(Mod(fft_out[2 : (N/2 + 1)]))

omeg <- 2*pi/L
omeg

model2 <- lm(y ~ sin(omeg*x) + cos(omeg*x) +sin(2*omeg*x)+ 
cos(2*omeg*x),data = DataC)
summary(model2)
x <- seq(DataC$x[1], DataC$x[N], by = 0.01)
myfit <- data.frame(x = x)
myfit$y <- predict(model2, newdata = myfit)
ggplot(DataC, aes(x, y)) +
geom_point() +
ggtitle(COUNTRY)+
geom_line(data = myfit, color = "blue")

###-----the coefficients for the Fourier series-----------------------
b11 = as.numeric(model2$coefficient[1]) #extracting the numeric valuess of the coefficients
b21 = as.numeric(model2$coefficients[2])
b31 = as.numeric(model2$coefficients[3])
b41 = as.numeric(model2$coefficients[4])
b51 = as.numeric(model2$coefficients[5])

roh=TestsP/(mean(prop_imm*Pop))*1 # the reporting rate



## measurement model C snippets

rmeas <-  Csnippet("
                   cases1 = rnbinom_mu((theta1),rho*H1);
                   cases2 = rnbinom_mu((theta2),rho*H2);")



dmeas <- Csnippet("
                 lik =  dnbinom_mu(cases1, (theta1), rho*H1, give_log) +
                  dnbinom_mu(cases2, (theta2), rho*H2, give_log);
                  ")


# define random simulator of measurement
rmeas.deaths <- " 
  deaths = rnbinom_mu((lam), CFR1 * D_new);
 "

# define evaluation of model prob density function
dmeas.deaths <- "
  lik = dnbinom_mu(deaths, (lam), CFR1 * D_new, give_log);
"

## initializer
rinit <- "
  double m = popsize/(S1_0 + E1_0 + I1_0 + R1_0 + D1_0+S2_0 + E2_0 + I2_0 + R2_0 + D2_0);
  S1 = nearbyint(m*S1_0);
  E1 = nearbyint(m*E1_0);
  I1 = nearbyint(m*I1_0);
  R1 = nearbyint(m*R1_0);
  D1 = nearbyint(m*D1_0);
  H1 = 0;
  S2 = nearbyint(m*S2_0);
  E2 = nearbyint(m*E2_0);
  I2 = nearbyint(m*I2_0);
  R2 = nearbyint(m*R2_0);
  D2 = nearbyint(m*D2_0);
  H2 = 0;
  D_new = 0;
  Phi = 0;
  noise = 0;
"

## rprocess
seas.sir.step2 <- "
  double rate[20];
  double dN[20];
  double Beta;
  double dW;
  double I;
  double Pop;
  I = I1 + I2;
  double immun;
  Pop = S1 + E1 + I1 + I2 + R1 + D1 + S2 + E2 + R2 + D2; //class 1 is the not immunized while 2 is the immunized
  Beta = exp(b1 + b2 * sin(M_2PI * Phi) + b3 * cos(M_2PI * Phi)+ b4 * sin(M_2PI * Phi*2) + b5 * cos(M_2PI * Phi*2));
  rate[0] = birthrate; // birth into susceptible1 class
  rate[1] =immunized; // rate of immunization
  rate[2] = contact_rate*severity_var*Beta*(I + iota) / (popsize); // force of infection1
  rate[3] = (1-infect_reduce)*contact_rate*severity_var*Beta * (I + iota) / (popsize); // force of infection2
  rate[4] = mu; // death from susceptible class1
  rate[5] = mu; // death from susceptible class2
  rate[6] = gamma1; // recovery rate1
  rate[7] = gamma2; // recovery rate2
  rate[8] = mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class1
  rate[9] = (1-immun/(popsize))*(1-mort_reduce)*mu2*Beta*(I + iota) / (popsize); // covid deaths from infectious class2
  rate[10] = mu; // death from recovered class1
  rate[11] = mu; // death from recovered class2
  rate[12] = mu; // death from I_1
  rate[13] = mu; // death from I_2
  rate[14] = lambda; // exposure (incubation period/latent period)1
  rate[15] = lambda; // exposure (incubation period/latent period)2
  rate[16] = epsilon;           // recovery to susceptible1
  rate[17] = epsilon*epsilonh_reduce;           // recovery to susceptible2
  rate[18] = mu;  //death from exposed1
  rate[19] = mu;  //death from exposed2

  dN[0] = rpois(rate[0] * dt); //births are Poisson
  reulermultinom(3, S1, &rate[14], dt, &dN[14]);
  reulermultinom(2, E1, &rate[2], dt, &dN[2]);
  reulermultinom(3, I1, &rate[6], dt, &dN[6]);
  reulermultinom(2, R1, &rate[16], dt, &dN[16]);
  reulermultinom(4, S2, &rate[15], dt, &dN[15]);
  reulermultinom(2, E2, &rate[3], dt, &dN[3]);
  reulermultinom(3, I2, &rate[7], dt, &dN[7]);
  reulermultinom(2, R2, &rate[17], dt, &dN[17]);

  // gamma noise (extra-demographic stochasticity)
  dW = rgammawn(sigmaSE,dt);

  S1 += dN[0] - dN[1] - dN[14] - dN[4];
  E1 += dN[14] - dN[2] - dN[18];
  I1 += dN[2] - dN[8] - dN[12]-dN[6];
  R1 += dN[6] - dN[10] - dN[16];
  D1 += dN[8];
  S2 += dN[1] + dN[16] + dN[17] - dN[15] - dN[5];
  E2 += dN[15] - dN[3] - dN[19];
  I2 += dN[3] - dN[7] - dN[9]-dN[13];
  R2 += dN[7] - dN[11] - dN[17];
  D2 += dN[9];
  Phi += dW;
  H1 += dN[2];
  H2 += dN[3];
  D_new += dN[8]+dN[9];
  noise += (dW - dt) / sigmaSE;
"
times=seq(0,2.83,by=1/52)
simulate(
  times=seq(0,2.83,by=1/52), t0=-1/52,
  covar=covariate_table(immun, times="time"),
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
   rinit=Csnippet(rinit),
  rprocess = euler(
    step.fun = Csnippet(seas.sir.step2), 
    delta.t = 1/52/20
  ),
  partrans=parameter_trans(
    toEst=Csnippet("
      to_log_barycentric(&T_S1_0,&S1_0,3);
      to_log_barycentric(&T_S2_0,&S2_0,3);
      T_sigmaSE = log(sigmaSE);
      T_iota = log(iota);"),
    fromEst=Csnippet("
      from_log_barycentric(&S1_0,&T_S1_0,3);
      from_log_barycentric(&S2_0,&T_S2_0,3);
      sigmaSE = exp(T_sigmaSE);
      iota = exp(T_iota);")
  ),
  statenames = c("S1","S2","E1","E2", "I1","I2", "R1","R2", "D1","D2","H1","H2", "Phi","noise","D_new"),
  obsnames=c("cases1","cases2"),
  paramnames = c("birthrate","gamma1","gamma2", "mu", "mu2","theta1", "theta2", "b1", "b2", "b3","b4","b5", "CFR",
    "popsize","rho", "iota", "sigmaSE","lambda","epsilon", "S1_0","E1_0", "I1_0", 
    "R1_0","D1_0","S2_0","E2_0", "I2_0","R2_0","D2_0","contact_rate",
"sev_reduce","severity_var","infect_reduce","epsilonh_reduce","mort_reduce"),
accumvars = c("H1","H2","D_new"),
params = c(popsize = Pop*(1-prop_pop_imm_21), birthrate=Anbirthrate,  lambda=lambdah,iota = iotah,
    gamma1 = gammah1,gamma2 = gammah2,
    mu = muh, mu2=CFR1, lam=lam, CFR=CFR1,b1 = b11, b2 = b21,b4=b41,b5=b51,
    b3 = b31, rho =roh, theta1 = thetah1, theta2 = thetah2,epsilon=epsilonhOmi, 
    sigmaSE = 0.03, S1_0 = 0.55 , E1_0 = 0.01 ,  I1_0 =0.01 , R1_0 = 0.01, D1_0 = 0.0,
    S2_0 = 0.4 , E2_0 = 0.01 ,  I2_0 =0.01 , R2_0 = 0.01, D2_0 = 0.0,
    contact_rate=contact_rateh,severity_var=severity_var,sev_reduce=sev_reduce,
    infect_reduce=infect_reduce,epsilonh_reduce=epsilonh_reduce,mort_reduce=mort_reduce),seed = 1914679908L) -> sir2C6AA

#plot(sir2C6AA)
sir2C6AAdta<-as(sir2C6AA,"data.frame")
head(sir2C6AAdta)

write.csv(sir2C6AAdta, "sir2C6AAdtaKE.csv")


##################################


times=seq(0,2.83,by=1/52)
tt=seq(1,length(times),by=1)

Cases = sir2C6AAdta$cases1 + sir2C6AAdta$cases2
S = sir2C6AAdta$S1 + sir2C6AAdta$S2
E = sir2C6AAdta$E1 + sir2C6AAdta$E2
I = sir2C6AAdta$I1 + sir2C6AAdta$I2
R = sir2C6AAdta$R1 + sir2C6AAdta$R2
D = sir2C6AAdta$D1 + sir2C6AAdta$D2
H = sir2C6AAdta$H1 + sir2C6AAdta$H2
DNs = sir2C6AAdta$D_new

Date=seq(as.Date("2022/1/1"), by = "week", length.out = length(times))
States=data.frame(COUNTRY,Date,Cases,S,E,I,R,D,H)
write.csv(States, "Current.csv")

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

DNH=DNs/deathreported
Dlow=DNs-1.96*sd(DNs)/sqrt(length(DNs))
Dlow[Dlow<0]<-0
Dhigh=DNs+1.96*sd(DNs)/sqrt(length(DNs))
DHlow=Dlow/deathreported
DHhigh=Dhigh/deathreported


lowH=H-1.96*sd(H)/sqrt(length(H))
highH=H+1.96*sd(H)/sqrt(length(H))

lowC=CasesCum-1.96*sd(CasesCum)/sqrt(length(CasesCum))
lowC[lowC<0]<-0
highC=CasesCum+1.96*sd(CasesCum)/sqrt(length(CasesCum))

lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))





Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


#Case=melt(CasesH, id.vars=c("Date"))
#head(Case)

###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())

low=Cases-1.96*sd(Cases)/sqrt(length(Cases))
low[low<0]<-0
high=Cases+1.96*sd(Cases)/sqrt(length(Cases))

lowH=H-1.96*sd(H)/sqrt(length(H))
lowH[lowH<0]<-0
highH=H+1.96*sd(H)/sqrt(length(H))


lowHC=Hcum-1.96*sd(Hcum)/sqrt(length(Hcum))
lowHC[lowHC<0]<-0
highHC=Hcum+1.96*sd(Hcum)/sqrt(length(Hcum))


Case=data.frame(Date,Cases,low,high)
Hidden=data.frame(Date,H,lowH,highH)

CaseC=data.frame(Date,CasesCum,lowC,highC)
HiddenC=data.frame(Date,Hcum,lowHC,highHC)


###ploting

library(lubridate)#for the dates
theme_set(theme_minimal())


C=ggplot(data = Case, aes(x = Date, y = CasesCum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Cases")+
geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = 0.1)
C
#png(file="C.png",width=600, height=350)

CC=ggplot(data = CaseC, aes(x = Date, y = Cases,group = 1))+
geom_line(color = "blue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
xlab("Time in Months") + ylab("Modelled Cases")+
geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1)
CC
#CCI=png(file="CC.png",width=600, height=350) 

HH=ggplot(data = HiddenC, aes(x = Date, y = Hcum,group = 1))+
geom_line(color = "red", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
xlab("Time in Months") + ylab("Cumulative Hidden Cases")+
geom_ribbon(aes(ymin = lowHC, ymax = highHC), alpha = 0.1)
HH
#write.xlsx(HH, "HHH.xlsx", append=TRUE,sheetName = COUNTRY)

Hp=ggplot(data = Hidden, aes(x = Date, y = H,group = 1))+
geom_line(color = "lightblue", size = 1)+
scale_x_date(date_breaks = "1 month",date_labels = "%b-\n%Y")+
geom_ma(ma_fun = SMA, n = 7, color = "red") + # Plot 7-day SMA
 xlab("Time in Months") + ylab("Hidden Cases")+
geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = 0.1)
Hp


sum(Cases) 
sum(H)

#####-------------Reporting-------------------

severity=c(0.88,0.08,0.03,0.01)
Sev=Cases%*%t(severity)

#library(data.table)
#library("Hmisc")
#library(xlsx)


colnames(Sev)<-c('Asymptomatic', 'Moderate', 'Severec', 'Critical')
Sev1=as.data.frame(Sev)
attach(Sev1)
lowAsy=Asymptomatic-1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))
lowAsy[lowAsy<0]<-0 
highAsy=Asymptomatic+1.96*sd(Asymptomatic)/sqrt(length(Asymptomatic))

lowMod=Moderate-1.96*sd(Moderate)/sqrt(length(Moderate))
lowMod[lowMod<0]<-0
highMod=Moderate+1.96*sd(Moderate)/sqrt(length(Moderate))

lowSev=Severec-1.96*sd(Severec)/sqrt(length(Severec))
lowSev[lowSev<0]<-0
higSev=Severec+1.96*sd(Severec)/sqrt(length(Severec))

lowCrit=Critical-1.96*sd(Critical)/sqrt(length(Critical))
lowCrit[lowCrit<0]<-0
highCrit=Critical+1.96*sd(Critical)/sqrt(length(Critical))


Sever=data.frame(Asymptomatic,lowAsy,highAsy,Moderate,lowMod, highMod,Severec,lowSev,higSev, Critical,lowCrit,highCrit,DNs,Dlow,Dhigh,DNH,DHlow,DHhigh)

StatesFull=data.frame(COUNTRY,Date,Cases,low,high,S,E,I,R,D,H,lowH,highH,Sever,lowHC, Hcum,highHC, lowC,CasesCum, highC)
head(StatesFull)

library(openxlsx) ##or
library(readxl)
library(WriteXLS)
write.xlsx2(StatesFull,"CurrentTrans200.xlsx", sheetName=COUNTRY, append=TRUE)

}


#transmissibility datasets


sheet = excel_sheets("CurrentTrans80.xlsx")
Combined_country_dataTrans80 = lapply(setNames(sheet, sheet), function(x) read_excel("CurrentTrans80.xlsx", sheet=x))
Combined_country_dataTrans80 = bind_rows(Combined_country_dataTrans80, .id="Sheet")
Combined_country_dataTrans80 =as.data.frame(Combined_country_dataTrans80)
Combined_country_dataTrans80$Date <- as.Date(Combined_country_dataTrans80$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_dataTrans80, "Combined_country_dataTrans80.xlsx")
head(Combined_country_dataTrans80)



sheet = excel_sheets("CurrentTrans120.xlsx")
Combined_country_dataTrans120 = lapply(setNames(sheet, sheet), function(x) read_excel("CurrentTrans120.xlsx", sheet=x))
Combined_country_dataTrans120 = bind_rows(Combined_country_dataTrans120, .id="Sheet")
Combined_country_dataTrans120 =as.data.frame(Combined_country_dataTrans120)
Combined_country_dataTrans120$Date <- as.Date(Combined_country_dataTrans120$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_dataTrans120, "Combined_country_dataTrans120.xlsx")
head(Combined_country_dataTrans120)


sheet = excel_sheets("CurrentTrans200.xlsx")
Combined_country_dataTrans200 = lapply(setNames(sheet, sheet), function(x) read_excel("CurrentTrans200.xlsx", sheet=x))
Combined_country_dataTrans200 = bind_rows(Combined_country_dataTrans200, .id="Sheet")
Combined_country_dataTrans200 =as.data.frame(Combined_country_dataTrans200)
Combined_country_dataTrans200$Date <- as.Date(Combined_country_dataTrans200$Date,format = "%m/%d/%Y")
write.xlsx (Combined_country_dataTrans200, "Combined_country_dataTrans200.xlsx")
head(Combined_country_dataTrans200)



