library(quantreg)
library(splines)

SEIR.model<-function (init, beta.s, gamma.e, gamma.i, times) {
  library(deSolve)
  seir <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- -beta.s * S * I
      dE <- beta.s * S * I - gamma.e * E
      dI <- gamma.e * E - gamma.i * I
      dR <- gamma.i * I
      return(list(c(dS, dE, dI, dR)))
    })
  }
  parameters <- c(beta.s = beta.s, gamma.e = gamma.e, gamma.i = gamma.i)
  out <- ode(y = init, times = times, func = seir, parms = parameters)
  out <- as.data.frame(out)
  out$time <- NULL
  return(out)
}

plot.func <- function (Tm, beta.s, gamma.e, gamma.i, p, Y, X) {
  init <- c(S = 1 - p, E = p/2, I = p/2, R = 0)
  mod <- SEIR.model(init, beta.s = beta.s, gamma.e = gamma.e, gamma.i, 1:Tm)
  matplot(x = 1:Tm, y = mod, type = "l", xlab = "Time", ylab = "Susceptible, Exposed, Infectious, and Recovered", main = "SEIR Model", lwd = 1, lty = 1, bty = "l", col = 1:4)
  xmax <- quantile(X, 0.975)
  matplot(x = 1:Tm, y = t(X), ylim = c(0, xmax), type = "p", xlab = "Time", ylab = "RNA levels", main = "", pch = 1, col = 1)
  Yv <- as.vector(Y)
  Xv <- as.vector(X)
  BX <- bs(Xv, df = 10)
  fit1 <- rq(Yv ~ BX, 0.05)
  fit2 <- rq(Yv ~ BX, 0.25)
  fit3 <- rq(Yv ~ BX, 0.5)
  fit4 <- rq(Yv ~ BX, 0.75)
  fit5 <- rq(Yv ~ BX, 0.95)
  plot(Xv, Yv, xlab = "RNA levels", ylab = "Active COVID-19 Cases", col = "gray", xlim = c(0, xmax))
  lines(sort(Xv), cbind(1, BX[order(Xv), ]) %*% coef(fit1), col = "red")
  lines(sort(Xv), cbind(1, BX[order(Xv), ]) %*% coef(fit2), col = "green")
  lines(sort(Xv), cbind(1, BX[order(Xv), ]) %*% coef(fit3), col = "blue")
  lines(sort(Xv), cbind(1, BX[order(Xv), ]) %*% coef(fit4), col = "green")
  lines(sort(Xv), cbind(1, BX[order(Xv), ]) %*% coef(fit5), col = "red")
  abline(0, 0)
  abline(500, 0)
  abline(1000, 0)
  abline(1500, 0)
  abline(2000, 0)
  abline(2500, 0)
  abline(3000, 0)
  abline(3500, 0)
  abline(4000, 0)
  abline(4500, 0)
  abline(5000, 0)
  abline(0, 0, lty = "dashed", col = "gray")
  abline(500, 0, lty = "dashed", col = "gray")
  abline(1000, 0, lty = "dashed", col = "gray")
  abline(1500, 0, lty = "dashed", col = "gray")
  abline(2000, 0, lty = "dashed", col = "gray")
  abline(2500, 0, lty = "dashed", col = "gray")
  abline(3000, 0, lty = "dashed", col = "gray")
  abline(3500, 0, lty = "dashed", col = "gray")
  abline(4000, 0, lty = "dashed", col = "gray")
  abline(4500, 0, lty = "dashed", col = "gray")
  abline(5000, 0, lty = "dashed", col = "gray")
}

Sewage.sim<-function (Tm, beta.s, gamma.e, gamma.i, p, N, mu.V.max, sd.V.max, mu.V.20, sd.V.20, T.V.max, Ts, Temp, mu.tau0, sd.tau0, mu.Q, sd.Q, G.mean, G.sd) {
  init <- c(S = 1 - p, E = p/2, I = p/2, R = 0)
  mod <- SEIR.model(init, beta.s = beta.s, gamma.e = gamma.e, gamma.i, 1:Tm)
  S <- mod[, 1]
  I <- mod[, 3]
  NC <- rpois(1:Tm, N * beta.s * S * I)
  Tr <- round(1/gamma.e + 1/gamma.i)
  NCA <- c(rep(0, (Tr - 1)), NC)
  CI <- rep(-99, Tm)
  for (i in 1:Tm) {
    CI[i] <- sum(NCA[i:(Tr + i)])
  }
  Vmat <- NULL
  for (i in 1:Tm) {
    if (NC[i] > 0) {
      Vmat <- rbind(Vmat, cbind(matrix(0, nrow = NC[i], ncol = (i - 1)), Viral.load.sim(NC[i], 1:30, mu.V.max, sd.V.max, mu.V.20, sd.V.20, T.V.max,  G.mean, G.sd), matrix(0, nrow = NC[i], ncol = (Tm - i + 1))))
    }
  }
  Vmat <- Vmat[, 1:Tm]
  TV <- apply(Vmat, 2, sum)
  tau0 <- rnorm(1, mu.tau0, sd.tau0)
  Q <- rnorm(1, mu.Q, sd.Q)
  tau.star <- tau0 * Q^(-(Temp - 20)/10)
  OVL <- TV * (1/2)^(Ts/tau.star)
  return(cbind(TV, CI, OVL))
}

Viral.load.sim <-function (n, t, mu.V.max, sd.V.max, mu.V.20, sd.V.20, T.V.max, G.mean, G.sd) {
  max <- rnorm(n, mu.V.max, sd.V.max)
  min <- rnorm(n, mu.V.20, sd.V.20)
  G <- rnorm(n, G.mean, G.sd)
  load <- NULL
  for (i in 1:n) {
    load.g <- 10^(t * max[i]/T.V.max) * (t <= T.V.max) + (t > T.V.max) * 10^(max[i] - (max[i] - min[i])/(20) * (t - T.V.max))
    load <- rbind(load, G[i] * load.g)
  }
  return(load)
}

MC.COVID19.wastewater<-function (Sim, Tm, beta.s, gamma.e, gamma.i, p, N, mu.V.max, sd.V.max, mu.V.20, sd.V.20, T.V.max, Ts, Temp, mu.tau0, sd.tau0, mu.Q, sd.Q, G.mean, G.sd){
  par(mfcol = c(3, 1))
  TV <- NULL
  CI <- NULL
  OVL <- NULL
  for (i in 1:Sim) {
    res <- Sewage.sim(Tm, beta.s, gamma.e, gamma.i, p, N, mu.V.max, sd.V.max, mu.V.20, sd.V.20, T.V.max, Ts, Temp, mu.tau0, sd.tau0, mu.Q, sd.Q, G.mean, G.sd)
    TV <- rbind(TV, res[, 1])
    CI <- rbind(CI, res[, 2])
    OVL <- rbind(OVL, res[, 3])
    print(i)
  }
  plot.func(Tm, beta.s, gamma.e, gamma.i, p, CI, OVL)
}
