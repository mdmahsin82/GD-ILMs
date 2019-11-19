###### SEIR data management and presentation of results ###
## Author: Md Mahsin
## Date: September 26, 2019


################################################

## Read the data ##

data <- read.table("SEIR_sim1.txt")
data <- read.table("SEIR_sim2.txt")
data <- read.table("SEIR_sim3.txt")
data <- read.table("SEIR_sim4.txt")
data <- read.table("SEIR_sim5.txt")
data <- read.table("SEIR_sim6.txt")
data <- read.table("SEIR_sim7.txt")
data <- read.table("SEIR_sim8.txt")
data <- read.table("SEIR_sim9.txt")
data <- read.table("SEIR_sim10.txt")

## Disctretize continuous infection time  based on minimum and maximum value 
## of the infection time

rw <- data[which.min(data$V4), ]
int_row <- rw$V1
##
data$V10 <- cut(data$V4, breaks= c(seq(0.9, 2.6, 0.1)), labels = c(2:18))

data$V10 <- cut(data$V4, breaks= c(seq(0.9, 2.71, 0.1)), labels = c(2:19))

data$V10 <- cut(data$V4, breaks= c(seq(0.9, 2.4, 0.1)), labels = c(2:16))

data$V10 <- cut(data$V4, breaks= c(seq(0.9, 2.9, 0.1)), labels = c(2:21))

data$V10 <- cut(data$V4, breaks= c(seq(0.9, 2.7, 0.1)), labels = c(2:19))

data$V10 <- cut(data$V4, breaks= c(seq(0.9, 2.62, 0.1)), labels = c(2:18))

data$V10 <- cut(data$V4, breaks= c(seq(0.98, 2.80, 0.1)), labels = c(2:19))

data$V10 <- cut(data$V4, breaks= c(seq(0.90, 2.80, 0.1)), labels = c(2:20))

data$V10 <- cut(data$V4, breaks= c(seq(0.90, 2.50, 0.1)), labels = c(2:17))

data$V10 <- cut(data$V4, breaks= c(seq(0.95, 2.55, 0.1)), labels = c(2:17))

data$V10 <- cut(data$V4, breaks= c(seq(0.90, 2.80, 0.1)), labels = c(2:20))
###
sdata <- subset(data, select = c(V1, V7, V8, V10))

tdata <- read.table("fclaimdata3.txt")

## Merge the data in R ##

mdata <- merge( tdata, sdata, by = "V1")

## make final data ##
shy <- subset(mdata, select = c(V1, V2, V3, V4, V10, V6))

shy$V10 <- as.numeric(as.character(shy$V10))

## Replace the minimum value as a initial infection ##

shy[int_row, 5] <- 1

write.table(shy, file = "seir_data1.txt", row.names = F, col.names = F)

##############

### After running the MCMC code in fortran, we used for suumarizing the results ##
######## Summarized MCMC outout for power-law kernel ####
## True parameter: alpha = 0.30, alpha1 = 0.30, beta = 2.0, lambda = 0.80, sigma = 0.60
## Summary of marginal posterior mean and 95% credible interval ##
## Fixed Effects posterior mcmc output ##
## Required packages ##
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library(doBy)
#######
myfun1 <- function(x){c(m=mean(x), v=var(x), lq = quantile(x, 0.025), uq = quantile(x, 0.975))}

## Data 1 ##
fixed <- read.table("parameter1.txt")
fixed1 <- fixed[-c(1:6000), - 6]
names(fixed1) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed1$iter <- 1:25000
ldat1 <- melt(fixed1, id.vars = "iter")
tt1 <- summaryBy(value ~ variable, data = ldat1, FUN =  myfun1)
tt1$data <- 1

## Data 2 ##
fixed <- read.table("parameter2.txt")
fixed2 <- fixed[-c(1:6000), -6]
names(fixed2) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed2$iter <- 1:25000
ldat2 <- melt(fixed2, id.vars = "iter")
tt2 <- summaryBy(value ~ variable, data = ldat2, FUN =  myfun1)
tt2$data <- 2

## Data 3 ##
fixed <- read.table("parameter3.txt")
fixed3 <- fixed[-c(1:6000), -6]
names(fixed3) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed3$iter <- 1:25000
ldat3 <- melt(fixed3, id.vars = "iter")
tt3 <- summaryBy(value ~ variable, data = ldat3, FUN =  myfun1)
tt3$data <- 3

## Data 4 ##
fixed <- read.table("parameter4.txt")
fixed4 <- fixed[-c(1:6000), -6]
names(fixed4) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed4$iter <- 1:25000
ldat4 <- melt(fixed4, id.vars = "iter")
tt4 <- summaryBy(value ~ variable, data = ldat4, FUN =  myfun1)
tt4$data <- 4

## Data 5 ##
fixed <- read.table("parameter5.txt")
fixed5 <- fixed[-c(1:6000), -6]
names(fixed5) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed5$iter <- 1:25000
ldat5 <- melt(fixed5, id.vars = "iter")
tt5 <- summaryBy(value ~ variable, data = ldat5, FUN =  myfun1)
tt5$data <- 5

## Data 6 ##
fixed <- read.table("parameter6.txt")
fixed6 <- fixed[-c(1:6000), -6]
names(fixed6) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed6$iter <- 1:25000
ldat6 <- melt(fixed6, id.vars = "iter")
tt6 <- summaryBy(value ~ variable, data = ldat6, FUN =  myfun1)
tt6$data <- 6

## Data 7 ###
fixed <- read.table("parameter7.txt")
fixed7 <- fixed[-c(1:6000), -6]
names(fixed7) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed7$iter <- 1:25000
ldat7 <- melt(fixed7, id.vars = "iter")
tt7 <- summaryBy(value ~ variable, data = ldat7, FUN =  myfun1)
tt7$data <- 7

## Data 8 ##
fixed <- read.table("parameter8.txt")
fixed8 <- fixed[-c(1:6000), -6]
names(fixed8) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed8$iter <- 1:25000
ldat8 <- melt(fixed8, id.vars = "iter")
tt8 <- summaryBy(value ~ variable, data = ldat8, FUN =  myfun1)
tt8$data <- 8

## Data 9 ##
fixed <- read.table("parameter9.txt")
fixed9 <- fixed[-c(1:6000), -6]
names(fixed9) <- c("alpha", "alpha1","delta", "lambda", "sigma")
fixed9$iter <- 1:25000
ldat9 <- melt(fixed9, id.vars = "iter")
tt9 <- summaryBy(value ~ variable, data = ldat9, FUN =  myfun1)
tt9$data <- 9

## Data 10 ##
fixed <- read.table("parameter10.txt")
fixed10 <- fixed[-c(1:6000), -6]
names(fixed10) <- c("alpha", "alpha1", "delta", "lambda", "sigma")
fixed10$iter <- 1:25000
ldat10 <- melt(fixed10, id.vars = "iter")
tt10 <- summaryBy(value ~ variable, data = ldat10, FUN =  myfun1)
tt10$data <- 10

###### Merge all the the summarized data ##
sumdata <- rbind(tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9, tt10)
names(sumdata) <- c("parameter", "pmean", "pvar", "plcd", "pucd", "data")
## Make it as a factor ##
sumdata$data <- as.factor(sumdata$data)

## For expression in the facet plot ##
levels(sumdata$parameter) <- expression("alpha", "alpha[1]", "delta" , "lambda", "sigma")
## To make true value line create a variable ##
sumdata$hline <- mapvalues(sumdata$parameter, from = c("alpha", "alpha[1]", "delta" , "lambda", "sigma"), 
to = c(0.3, 0.3, 2.0, 0.85, 0.90))

## Make it as continuous ##
sumdata$hline <- as.numeric(as.character(sumdata$hline))

## Make a plot for all the data ##
pd <- position_dodge(0.50)

splot <- ggplot(sumdata, aes(x = data, y = pmean, group = data)) + geom_point(size = 0.4, position=pd) + 
geom_errorbar(data = sumdata, aes(ymin = plcd, ymax = pucd), width=0.2, size = 0.2, position=pd) + 
facet_wrap(~ parameter, scales="free", labeller = label_parsed, ncol=2)  + geom_hline(aes(yintercept = hline), 
size = 0.2, color = "red") + xlab("") + ylab("") + theme_bw() + theme(axis.title.x=element_blank(), 
axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="none", panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(), strip.text.x = element_text(margin = margin(0.8, 0.8, 0.8, 0.8)))

ggsave(splot, file = "fixed_seir.pdf")

###### Random effects ##

#phi <- scan("phi1.txt")
#phi <- scan("phi2.txt")
#phi <- scan("phi3.txt")
#phi <- scan("phi4.txt")
#phi <- scan("phi5.txt")
#phi <- scan("phi6.txt")
#phi <- scan("phi7.txt")
#phi <- scan("phi8.txt")
#phi <- scan("phi9.txt")
phi <- scan("phi10.txt")
nplot<- 31000
wdim <- 16
phi.mtx <- matrix(0,nplot,wdim)
for (i in 1:nplot) {
    for (j in 1:wdim) {
       phi.mtx[i,j]=phi[(i-1)*wdim+j]
    }
 }
##########

### Data 1 ##
phi1 <- phi.mtx[-c(1:6000), ]
phi1 <- as.data.frame(phi1)
names(phi1) <- paste("phi", 1:16, sep = "")
phi1$iter <- 1:25000
lphi1 <- melt(phi1, id.vars = "iter")
tt1 <- summaryBy(value ~ variable, data = lphi1, FUN =  myfun1)
tt1$data <- 1

#### Data 2 ##

phi2 <- phi.mtx[-c(1:6000), ]
phi2 <- as.data.frame(phi2)
names(phi2) <- paste("phi", 1:16, sep = "")
phi2$iter <- 1:25000
lphi2 <- melt(phi2, id.vars = "iter")
tt2 <- summaryBy(value ~ variable, data = lphi2, FUN =  myfun1)
tt2$data <- 2

#### Data 3 ##

phi3 <- phi.mtx[-c(1:6000), ]
phi3 <- as.data.frame(phi3)
names(phi3) <- paste("phi", 1:16, sep = "")
phi3$iter <- 1:25000
lphi3 <- melt(phi3, id.vars = "iter")
tt3 <- summaryBy(value ~ variable, data = lphi3, FUN =  myfun1)
tt3$data <- 3


#### Data 4 ##

phi4 <- phi.mtx[-c(1:6000), ]
phi4 <- as.data.frame(phi4)
names(phi4) <- paste("phi", 1:16, sep = "")
phi4$iter <- 1:25000
lphi4 <- melt(phi4, id.vars = "iter")
tt4 <- summaryBy(value ~ variable, data = lphi4, FUN =  myfun1)
tt4$data <- 4

#### Data 5 ##

phi5 <- phi.mtx[-c(1:6000), ]
phi5 <- as.data.frame(phi5)
names(phi5) <- paste("phi", 1:16, sep = "")
phi5$iter <- 1:25000
lphi5 <- melt(phi5, id.vars = "iter")
tt5 <- summaryBy(value ~ variable, data = lphi5, FUN =  myfun1)
tt5$data <- 5

#### Data 6 ##

phi6 <- phi.mtx[-c(1:6000), ]
phi6 <- as.data.frame(phi6)
names(phi6) <- paste("phi", 1:16, sep = "")
phi6$iter <- 1:25000
lphi6 <- melt(phi6, id.vars = "iter")
tt6 <- summaryBy(value ~ variable, data = lphi6, FUN =  myfun1)
tt6$data <- 6

#### Data 7 ##

phi7 <- phi.mtx[-c(1:6000), ]
phi7 <- as.data.frame(phi7)
names(phi7) <- paste("phi", 1:16, sep = "")
phi7$iter <- 1:25000
lphi7 <- melt(phi7, id.vars = "iter")
tt7 <- summaryBy(value ~ variable, data = lphi7, FUN =  myfun1)
tt7$data <- 7

#### Data 8 ##

phi8 <- phi.mtx[-c(1:6000), ]
phi8 <- as.data.frame(phi8)
names(phi8) <- paste("phi", 1:16, sep = "")
phi8$iter <- 1:25000
lphi8 <- melt(phi8, id.vars = "iter")
tt8 <- summaryBy(value ~ variable, data = lphi8, FUN =  myfun1)
tt8$data <- 8

#### Data 9 ##

phi9 <- phi.mtx[-c(1:6000), ]
phi9 <- as.data.frame(phi9)
names(phi9) <- paste("phi", 1:16, sep = "")
phi9$iter <- 1:25000
lphi9 <- melt(phi9, id.vars = "iter")
tt9 <- summaryBy(value ~ variable, data = lphi9, FUN =  myfun1)
tt9$data <- 9

#### Data 10 ##

phi10 <- phi.mtx[-c(1:6000), ]
phi10 <- as.data.frame(phi10)
names(phi10) <- paste("phi", 1:16, sep = "")
phi10$iter <- 1:25000
lphi10 <- melt(phi10, id.vars = "iter")
tt10 <- summaryBy(value ~ variable, data = lphi10, FUN =  myfun1)
tt10$data <- 10


###### Merge all the the summarized data ##
sumdata <- rbind(tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9, tt10)
names(sumdata) <- c("parameter", "pmean", "pvar", "plcd", "pucd", "data")

ssdata <- subset(sumdata, parameter == "phi1" | parameter == "phi2" |parameter == "phi3" |
			parameter == "phi4" |parameter == "phi5" |parameter == "phi6")
## Make it as a factor ##
sumdata$data <- as.factor(sumdata$data)

ssdata$data <- as.factor(ssdata$data)

## For expression in the facet plot ##


levels(ssdata$parameter) <- expression("phi[1]", "phi[2]", "phi[3]", "phi[4]", 
					"phi[5]", "phi[6]", "phi[7]", "phi[8]",
					"phi[9]", "phi[10]", "phi[11]", "phi[12]",
					"phi[13]", "phi[14]", "phi[15]", "phi[16]")

# To make true value line create a variable ##

ssdata$hline <- mapvalues(ssdata$parameter, from = c("phi[1]", "phi[2]", "phi[3]", "phi[4]",
				"phi[5]", "phi[6]"),
				to = c(-0.317, -0.011, -0.396, -0.249, -0.632,-1.021))

## Make it as continuous ##
sumdata$hline <- as.numeric(as.character(sumdata$hline))

ssdata$hline <- as.numeric(as.character(ssdata$hline))

## Make a plot for all the data ##
pd <- position_dodge(0.50)


splot3 <- ggplot(ssdata, aes(x = data, y = pmean, group = data)) + geom_point(size = 0.4, position=pd) + 
geom_errorbar(data = ssdata, aes(ymin = plcd, ymax = pucd), width=0.2, size = 0.2, position=pd) + 
facet_wrap(~ parameter,labeller = label_parsed, ncol=2)  + geom_hline(aes(yintercept = hline), size = 0.2, color = "red") + 
xlab("") + ylab("") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
axis.ticks.x=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
strip.text.x = element_text(margin = margin(0, 0, 0, 0)))




ggsave(splot3, file = "random_seir6.pdf", scale = 1.1)

#######################




