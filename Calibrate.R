
########################################################################
#                          Calibrate function                          #
########################################################################


# Author: Ali Rafei (arafei {AT} umich {DOT} edu)
#
# calibrate(psmp, nsmp, method=c("PSDO", "PMM"), model=c("GLM1", "GLM2", "BART", "CART", "RF"), q=0.01, JRR=FALSE, JKn=c(100, 100))
# Generates a set of pseudo-weights for the non-probability sample and a set of replicate pseudo-weights based on a modified Jackknife method
# Arguments:
# psmp: Probability sample with auxiliary variables and weights in the last column
# nsmp: Nonprobability sample with auxiliary variables
# method: calibrating method, either PSDO or PMM
# model: The model should be used for calibration: GLM, BART, CART or Random Forest
# q: The number of subclasses in PMM: 0.01 means to use percentiles of the distribution
# JRR: logical value to determine if Jackknife replicate method is required
# JKn0: favorite number of replicates for the probability sample
# JKn1: favorite number of replicates for the nonprobability sample

# Required packages:
library(MASS)
library(pps)
library(BayesTree)
require(BART, quietly=TRUE)
library(survey)
require(rpart)
require(randomForest, quietly=TRUE)

calibrate <- function(psmp, nsmp, method=c("PSDO", "PMM"), model=c("GLM1", "GLM2", "BART", "CART", "RF"), q=0.01, JRR=FALSE, JKn=c(100, 100)){
	out <- calsmp(psmp1=psmp, nsmp1=nsmp, method1=method, model1=model, q1=q)
	if(JRR){
		n0 <- nrow(psmp)
		n1 <- nrow(nsmp)
		out <- data.frame(psdo_wght=out, matrix(0, n1, JKn[1]+JKn[2]))
		names(out)[-1] <- paste("psdo", 1:(JKn[1]+JKn[2]), sep="")
		jk_id <- c(rep(1:JKn[1], rep(n0%/%JKn[1], JKn[1])), rep(JKn[1], n0%%JKn[1]))
		for(j in 1:JKn[1]){
			tmp <- calsmp(psmp1=psmp[jk_id!=j, ], nsmp1=nsmp, method=method1, model1=model, q=q1)
			out[, paste("psdo", j, sep="")] <- ((n0-sum(jk_id==j))/n0)*tmp
		}
		jk_id <- c(rep((JKn[1]+1):(JKn[1]+JKn[2]), rep(n1%/%JKn[2], JKn[2])), rep((JKn[1]+JKn[2]), n1%%JKn[2]))	
		for(j in (JKn[1]+1):(JKn[1]+JKn[2])){
			tmp <- calsmp(psmp1=psmp, nsmp1=nsmp[jk_id!=j, ], method1=method, model1=model, q=q)
			out[jk_id!=j, paste("psdo", j, sep="")] <- ((n0-sum(jk_id==j))/n0)*tmp		
		}
	}
	out
}

calsmp <- function(psmp1, nsmp1, method1=c("PSDO", "PMM"), model1=c("GLM1", "GLM2", "BART", "CART", "RF"), q1=0.01){
	pvar <- ncol(nsmp1)
	smp <- rbind(psmp1[, 1:pvar], nsmp1)
	smp$z <- rep(0:1, c(nrow(psmp1), nrow(nsmp1)))
	psmp1$Wi_log <- log(1/psmp1[, pvar+1])
	if(model1=="GLM1"){
		fit1 <- lm(Wi_log~., psmp1[, -(pvar+1)])
		nsmp1$Wt <- predict(fit1, newdata=nsmp1)
		if(method1=="PSDO"){
			fit2 <- glm(z~., smp, family=binomial(link="logit"))
			nsmp1$PZi <- predict(fit2, newdata=nsmp1[, 1:pvar], type="response")
			wght <- exp(-nsmp1$Wt+log(1-nsmp1$PZi)-log(nsmp1$PZi))
		}
	}else if(model1=="GLM2"){
		fit1 <- lm(Wi_log~.^2, psmp1[, -(pvar+1)])
		nsmp1$Wt <- predict(fit1, newdata=nsmp1)
		if(method1=="PSDO"){
			fit2 <- glm(z~.^2, smp, family=binomial(link="logit"))
			nsmp1$PZi <- predict(fit2, newdata=nsmp1[, 1:pvar], type="response")
			wght <- exp(-nsmp1$Wt+log(1-nsmp1$PZi)-log(nsmp1$PZi))
		}
	}else if(model1=="BART"){
		X_smp <- makeind(smp[, 1:pvar])
		fit1 <- mc.wbart(x.train=X_smp[smp$z==0, ], y.train=psmp1$Wi_log, x.test=X_smp[smp$z==1, ], ntree=200L, mc.cores=192L, nskip=100L, ndpost=1000L)
		nsmp1$Wt <- fit1$yhat.test.mean
		if(method1=="PSDO"){
			fit2 <- mc.pbart(x.train=X_smp, y.train=smp$z, x.test=X_smp[smp$z==1, ], ntree=200L, mc.cores=192L, nskip=100L, ndpost=1000L)
			nsmp1$PZi <- apply(pnorm(fit2$yhat.test), 2, median)
			wght <- exp(-nsmp1$Wt+log(1-nsmp1$PZi)-log(nsmp1$PZi))
		}
	} else if(model1=="CART"){
		fit1 <- rpart(Wi_log~., data=psmp1[, -(pvar+1)], control=rpart.control(minbucket=50, cp=0.0,  xval=0), method="anova")
		nsmp1$Wt <- predict(fit1, newdata=nsmp1, type="vector")
		if(method1=="PSDO"){
			fit2 <- rpart(factor(z)~., smp, control=rpart.control(minbucket=50, cp=0.0, xval=0), method="class")
			nsmp1$PZi <- as.numeric(predict(fit2, newdata=nsmp1[, 1:pvar], type="prob")[, 2])
			wght <- exp(-nsmp1$Wt+log(1-nsmp1$PZi)-log(nsmp1$PZi))
		}
	} else if(model1=="RF"){
		fit1 <- randomForest(Wi_log~., data=psmp1[, -(pvar+1)], ntree=100, mtry=ceiling(sqrt(pvar)), importance=FALSE, na.action=na.omit, replace=FALSE)
		nsmp1$Wt <- predict(fit1, newdata=nsmp1, type="response")
		if(method1=="PSDO"){
			fit2 <- randomForest(factor(z)~., smp, ntree=100, mtry=ceiling(sqrt(pvar)), importance=FALSE, na.action=na.omit, replace=FALSE)
			nsmp1$PZi <- as.numeric(predict(fit2, newdata=nsmp1[, 1:pvar], type="prob")[, 2])
			wght <- exp(-nsmp1$Wt+log(1-nsmp1$PZi)-log(nsmp1$PZi))
		}
	}
	
	if(method1=="PMM"){
		p <- quantile(psmp1$Wi_log, seq(0, 1, by=q1))
		psmp1$Wi_q <- as.numeric(cut(psmp1$Wi_log, breaks=unique(p), levels=1:(length(unique(p))-1), labels=1:(length(unique(p))-1), include.lowest=T))
		nsmp1$Wi_q <- as.numeric(cut(nsmp1$Wt, breaks=unique(p), levels=1:(length(unique(p))-1), labels=1:(length(unique(p))-1), include.lowest=T))
		qk <- by(psmp1[, pvar+1], psmp1$Wi_q, function(x)sum(x, rm.na=T))
		tmp_qk <- data.frame(Wi_q=as.numeric(names(qk)), qk=as.numeric(qk))
		nsmp1$ID1 <- 1:nrow(nsmp1)
		nsmp1 <- merge(nsmp1, tmp_qk, by="Wi_q", all.x=T, sort=F)
		nsmp1 <- nsmp1[order(nsmp1$ID1), ]
		nsmp1$nk <- ave(nsmp1$Wi_q, nsmp1$Wi_q, FUN=length)
		nsmp1$nk_n <- ave(nsmp1$Wi_q, nsmp1$Wi_q, FUN=seq_along)
		PMM_qk <- nsmp1[nsmp1$nk_n==1, c("qk", "nk")]
		qk_sum <- sum(PMM_qk$qk, na.rm=T)
		nk_sum <- sum(PMM_qk$nk, na.rm=T)
		wght <- (nsmp1$qk/qk_sum)/(nsmp1$nk/nk_sum)
		nsmp1[, c("Wt", "PZi", "ID1", "Wi_q", "qk", "nk", "nk_n")] <- NULL
	}
	wght
}

