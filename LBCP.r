
library(nmcdr)
library(changepoint)
library(StepSignalMargiLike)
#####one change point
nsim <- 500
n <- 100
sigma <- 0.5
u0 <- 1
u1 <- 3
prifuns <- prifun2
prifuns1 <- prifun2

fix  <- 0
                                        #gamma <- rep(NA, n)
kp <- 2
tau <- 2
nu <- 1
pu <- 0
psig <- 10

prifun0 <- function(u, pprm){
    pu <- pprm[1]

    psig <- pprm[2]
    dnorm(u, pu, psig)
}
prifun1 <- function(u, pprm){
    kp <- pprm[1]
    pu <- pprm[2]
    psig <- pprm[3]
    f <- function(u)
        { (u)^(2* kp) * dnorm(u, pu, psig) }
                                        #  tauk <- integrate(f, -Inf, Inf, rel.tol = 1e-5, subdivisions = 1000L)$value
    (u)^{2*kp} *  dnorm(u, pu, psig)
}
prifun2 <- function(u, pprm){
    kp <- pprm[1]
    tau <- pprm[2]
    nu <- pprm[3]
    kp * tau^(nu/2)/gamma(nu/(2* kp))  * ((u)^2)^(-(nu+1)/2) * exp(-((u)^2/tau)^(-kp))
}




LBCP <- function(data, prior, mindis, pprm, kmax, mkn = NULL){
    if(prior == 'nonlcim'){
        prifuns <- prifun2
        kp <- pprm[1]
        tau <- pprm[2]
        nu <- pprm[3]
    }else if(prior == 'nonlcm'){
        prifuns <- prifun1
        kp <- pprm[1]
        pu <- pprm[2]
        psig <- pprm[3]
    }else{
        prifuns <- prifun0
        pu <- pprm[1]
        psig <- pprm[2]
    }
    Y <- (data - mean(data))/sd(data)
    n <- length(Y)
    mv <- rep(1, n)
    nI <- mindis
    ests <-  sqrt(1/2)
    integrand0 <- function(u, X21, prifun, pprm){
        lu <- length(u)
        mY <- matrix(rep(X21, lu), ncol = lu)
        mu <- matrix(rep(u, nrow(mY)), ncol = lu, byrow = T)
        apply((2 * pi * ests^2)^(-1/2) * exp(-(mY - mu)^2/ (2 *  ests^2)), 2, prod ) * prifun(u, pprm)
    }
    tstat <- rep(0, n)
    pmk <- rep(0, n)
    for(i in nI : (n - nI)){
        X1 <- Y[max((i - nI ), 1) : (i-1 )] 
        X2 <- Y[(i) :  (i + nI-1 )]
        X21 <- (X2 - mean(X1))
        tstat[i] <- (integrate(integrand0, -Inf, Inf, X21, prifuns, pprm)$value/prod(dnorm(X21, 0, ests)))
    }
    pmk <- tstat
    ot <- order(pmk, decreasing = T)
    ordtstat <- pmk[ot]
    candscreen <- (1 : n)[ot]
    newcand <- candscreen
    maxord <- ordtstat
    i <- 1
    resp <- 0
    while(length(resp) <= kmax & i <= n -nI){
        if(min(abs(newcand[i] - resp)) >= nI){
            resp <- c(resp, newcand[i])
        }
        i <- i+ 1
    }
      pmk <- tstat
            resp <- 0
            for(i in nI:(n- nI)){
                if(pmk[i] == max(pmk[(i - nI + 1) : (i + nI)]))
                    resp <- c(resp, i)
            }
    screenI <- c(1, resp[-1][order(resp[-1])])
    cadpt <- screenI
    sdY <-Ybar <- Y
    l <- length(cadpt)
    if(l > 2){
        for(i in 2 : (l - 1)){
            Ybar[cadpt[i] : (cadpt[i + 1] -1)]<- mean(Y[cadpt[i -1] : (cadpt[i] -1)])
        }
    }else{
        Ybar[cadpt[2] : n]<- mean(Y[cadpt[2 -1] : (cadpt[2] -1)])
    }
    Ybar[cadpt[l] : n] <- mean(Y[cadpt[l -1] : (cadpt[l] -1)])
    newY <- (Y- Ybar)
    newY[cadpt[1] : (cadpt[2] - 1)] <- (Y[cadpt[1] : (cadpt[2] - 1)] - mean(Y[cadpt[1] : (cadpt[2] - 1)]))
    integrand <- function(u, i, prifun, pprm){
        lu <- length(u)
        if(i + 1 <= l) umY <- cadpt[i + 1] -1 else umY <- length(newY)
        mY <- matrix(rep(newY[cadpt[i] : umY], lu), ncol = lu)
        mu <- matrix(rep(u, nrow(mY)), ncol = lu, byrow = T)
        apply((2 * pi * ests^2)^(-1/2) * exp(-(mY - mu)^2/ (2 *  ests^2)), 2, prod ) * prifun(u, pprm)
    }
    chpint <- res <-rep(NA, l)
    if(l > 2){
        for(i in 2:(l -1)){
            chpint[i] <- log(integrate(integrand, -Inf, Inf, i, prifuns, pprm)$value)
            res[i] <- sum(log(dnorm(newY[-c( (cadpt[i] : (cadpt[i + 1] -1)))], 0, ests)))  + chpint[i]
        }
        i <- l
        chpint[i] <- log(integrate(integrand, -Inf, Inf, i, prifuns, pprm)$value)
        res[i] <- sum(log(dnorm(newY[-(cadpt[i] : n)], 0, ests)))  + chpint[i]
    }else{
        chpint[l] <- log(integrate(integrand, -Inf, Inf, l, prifuns, pprm)$value)
        res[l] <- sum(log(dnorm(newY[-c( (cadpt[2] : (n)))], 0, ests)))  + chpint[l]
    }
    o <- order(res[-1], decreasing = T)
    flg <- 0
    if(!is.null(mkn)){
        selectp <- cadpt[-1][o][1 : mkn][order(cadpt[-1][o][1 : mkn])]
        return(list(lbcp = selectp, nbcp = mkn))
    }else{
        current <-  sum(log(dnorm(newY, 0, ests)))
        for(kn in 1 : (l - 1)){
            nonordered  <- cadpt[-1][o][1:kn]
            chpint[chpint == -Inf] <- -1e6
            tchpt <- c(1, nonordered, n+ 1)
            lt <- kn + 1
            dnormjump <- mu <- vchpint <- rep(NA, lt)
            jumpindex <- 0
            if(lt > 2){
                for(j in 2 :  (lt)){
                    ixcad <- which(cadpt == tchpt[j] )
                    vchpint[j] <- (chpint[ixcad])
                    jumpindex <- c(jumpindex, (tchpt[j] : min((cadpt[ixcad + 1] - 1), n, na.rm = T)))
                }
            }else{
                ixcad <- which(cadpt == tchpt[2] )
                vchpint[2] <- sum(chpint[ixcad])
                if(ixcad + 1 < l ){
                    jumpindex <- c(jumpindex, (tchpt[2] : (cadpt[ixcad + 1] - 1)))
                }else{
                    jumpindex <- c(jumpindex, (tchpt[2] : n))
                }
            }
            jumpindex <- jumpindex[-1]
            current <- c(current, sum((vchpint), na.rm = T) +  sum(log(dnorm(newY[-c( jumpindex)], 0, ests))))
        }
        selectp <- cadpt[-1][o]
        ix <- current == max(current)
        mkn <- min(which(ix)) - 1
        selectp <- selectp[1:mkn]
        selectp <- selectp[order(selectp)]
        return(list(lbcp = selectp, nbcp = mkn))
    }

}

predict.BCP <- function(Y, cpt) {
 
    
    cpt <- sort(unique(c(cpt, 0, length(Y))))
    
    
   
      fit <- rep(0, length(Y))
      
      for (i in 1:(length(cpt) - 1)) {
        fit[(cpt[i] + 1):cpt[i + 1]] <- mean(Y[(cpt[i] + 1):cpt[i + 1]])
        
      }
      
   
    
    return(fit)
    
}

plot.BCP <- function(Y,cpt, ...){
    plot(Y, type="l", ...)
    lines(x=predict.BCP(Y,cpt), type="l",col="red")  
    title("Data and the fitted signal")
}

plot.BCPcompare <- function(Y,cpt,cpt1, ylab, xaxt){
    plot(Y, type="p", ylab = ylab, xaxt = xaxt, lwd = 0.1, cex = 0.1)
    lines(x=predict.BCP(Y,cpt), type="l",col="red", lty = 1, lwd =2)
    lines(x=predict.BCP(Y,cpt1), type="l",col="blue", lty = 3, lwd = 2) 
   title(ylab = ylab, outer = TRUE, line = 3)
}

plot.BCP3 <- function(Y,cpt,cpt1, cpt2, ylab, xaxt){
    plot(Y, type="p", ylab = ylab, xaxt = xaxt, lwd = 0.1, cex = 0.1)
    lines(x=predict.BCP(Y,cpt), type="l",col="red", lty = 1, lwd =2)
    lines(x=predict.BCP(Y,cpt1), type="l",col="blue", lty = 6, lwd = 2)
    lines(x=predict.BCP(Y,cpt2), type="l",col="gold", lty = 3, lwd = 2) 
   
}
postscript('temp.ps')
par(mfrow = c(2, 1),  mar=c(2, 4, 1, 1) + 0.1)
plot.BCPcompare(aggratio[, 2], hospLBCP, hospnot, "Hospitalization rate", 'n')
plot.BCPcompare(mtemp[, 2], tempLBCP, tempnot, "Temperature", NULL)
dev.off()

plot.BCP3(aggratio[, 2], hospcpt, hospwbs, hospstep, "Hospitalization rate", 'n')
plot.BCP3(mtemp[, 2], tempcpt, tempwbs, tempstep, "Temperature", NULL)
#plot.BCP3(mtemp[, 2], tempLBCP, tempnot, "Temperature", NULL)


plot.BCPcompare(aggratio[, 2], hospLBCP, hospnot)
 plot.BCPcompare(aggratio[, 2], hospLBCP, hospnot)
hospLBCP <- LBCP(aggratio[, 2], 'nonlocim', 14, c(10, 2, 1), 50)[[1]]
tempLBCP <- (LBCP(mtemp[, 2], 'nonlocim', 14, c(10, 2, 1), 50)[[1]])
plot.BCP(aggratio[, 2], hospLBCP, ylab = 'Hospitalization rate')
plot.BCP(mtemp[, 2], tempLBCP, ylab = 'Temperature')

qqplot(hospLBCP, tempLBCP,  ylab= 'Temperature', xlab = 'Standardized log hospitalization rate' ,col = 2, pch = 19, cex.lab = 1.4)
abline(0, 1)
legend("topleft", "BCP", bty = "n", cex = 2)
 prior.B <- prior.norm.B(aggratio[, 2 ])
hospstep <- est.changepoints(aggratio[, 2 ], model = 'normal', prior.B, 100)

prior.B <- prior.norm.B(mtemp[, 2 ])
tempstep <- est.changepoints(mtemp[, 2 ], model = 'normal', prior.B, 100)
qqplot(hospstep, tempstep,  ylab= 'Temperature', xlab = 'Standardized log hospitalization rate ', col = 2, pch = 19, cex.lab = 1.4)
abline(0, 1) ##temp 52, hosp 68
legend("topleft", "SML", bty = "n", cex = 2)
tempcpt <- cpts(cpt.mean(mtemp[, 2],penalty="MBIC",pen.value=0,method="PELT",Q=kmax,test.stat="Normal",class=TRUE,  param.estimates=TRUE,minseglen=21))

hospcpt <- cpts(cpt.mean(aggratio[, 2 ],penalty="MBIC",pen.value=0,method="PELT",Q=100,test.stat="Normal",class=TRUE,  param.estimates=TRUE,minseglen=21))
qqplot(hospcpt, tempcpt,  ylab= 'Temperature', xlab = 'Standardized log hospitalization rate', col = 2, pch = 19, cex.lab = 1.4)
abline(0, 1)
legend("topleft", "PELT", bty = "n", cex = 2)
qqplot(hospnot, tempnot,  ylab= 'Temperature', xlab = 'Standardized log hospitalization rate ', col = 2, pch = 19, cex.lab = 1.4)
abline(0, 1)
legend("topleft", "NOT", bty = "n", cex = 2)
qqplot(hospwbs, tempwbs,  ylab= 'Temperature', xlab = 'Standardized log hospitalization rate ', col = 2, pch = 19 ,cex.lab = 1.4)
abline(0, 1)
legend("topleft", "WBS", bty = "n", cex = 2)
plotchange<- function(data, chgpt,  xlab = NULL, ylab = NULL, title = NULL, plot){
    Z<- data
    Z = data.frame(timestamp = 1:length(Z), count = Z)
    count <- x <- y <- xend <- yend <- 0L
    g = ggplot2::ggplot(Z, ggplot2::aes(x = timestamp, y = count)) + 
        ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), 
                                             panel.grid.major = ggplot2::element_blank())
    if (is.null(xlab)) 
        xlab = " "
    if (is.null(ylab)) 
        ylab = ""
    if (is.null(title)) 
        title = " "
    g = g + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)
    g = g + ggplot2::geom_line()

    v =  chgpt #temp[order(temp)]
    v = c(0, v)
    for (j in 2:length(v)) {
        M = mean(Z$count[(v[j - 1] + 1):v[j]])
        df2 = data.frame(Z$timestamp[v[j]], Z$timestamp[v[j]], -Inf, M)
        names(df2) = c("x", "xend", "y", "yend")
        g = g + ggplot2::geom_segment(data = df2, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linetype = 2, size = 0.5)
        g = g + ggplot2::guides(color = TRUE) 
    }
    
    
    
    g = g + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,  0)) + theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
    if(plot == TRUE){
        plot(g)
    }

}


hospnot <- features(not((aggratio[, 2] - mean(aggratio[, 2]))/sd(aggratio[, 2]),  M = 5000, method = c("not"),contrast = c("pcwsConstMean"), rand.intervals = nI,parallel = FALSE, augmented = FALSE, intervals))$cpt
hospwbs <- changepoints(wbs((aggratio[, 2] - mean(aggratio[, 2]))/sd(aggratio[, 2]), M = 5000))$cpt.th[[1]]


tempnot <- features(not((mtemp[, 2] - mean(mtemp[,2]))/sd(mtemp[, 2]),  M = 50000, method = c("not"),contrast = c("pcwsConstMean"), rand.intervals = TRUE,parallel = FALSE, augmented = FALSE, intervals))$cpt
tempwbs <- changepoints(wbs((mtemp[, 2] - mean(mtemp[,2]))/sd(mtemp[, 2]), M = 5000))$cpt.th[[1]]

pdf('hospchgp.pdf')
plot.BCP(aggratio[, 2], hospLBCP, NULL, ylab = 'Hospitalization rate')
dev.off()
pdf('tempchgp.pdf')
plot.BCP(mtemp[, 2], tempLBCP, NULL, ylab = 'Temperature')
dev.off()

pdf('hospcpt.pdf')
plot.BCP(aggratio[, 2], hospcpt, NULL, ylab = 'Hospitalization rate')
dev.off()
pdf('tempcpt.pdf')
plot.BCP(mtemp[, 2], tempcpt, NULL, ylab = 'Temperature')
dev.off()

pdf('stephospchgp.pdf')
plot.BCP(aggratio[, 2], hospstep, NULL, ylab = 'Hospitalization rate')
dev.off()
pdf('steptempchgp.pdf')
plot.BCP(mtemp[, 2], tempstep, NULL, ylab = 'Temperature')
dev.off()


pdf('nothospchgp.pdf')
plot.BCP(aggratio[, 2], hospnot, NULL, ylab = 'Hospitalization rate')
dev.off()
pdf('nottempchgp.pdf')
plot.BCP(mtemp[, 2], tempnot, NULL, ylab = 'Temperature')
dev.off()

pdf('wbshospchgp.pdf')
plot.BCP(aggratio[, 2], hospwbs, NULL, ylab = 'Hospitalization rate')
dev.off()
pdf('wbstempchgp.pdf')
plot.BCP(mtemp[, 2], tempwbs, NULL, ylab = 'Temperature')
dev.off()
