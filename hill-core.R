library(polynom)
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rgl)
library(plotrix)
library(calibrate)

# Hill
hill <- function(n, a, g, b = 1) {
  f <- function(x) a / (1 + x ^ g)
  df <- function(x) -a * g * x ^ (g - 1) / (x ^ g  + 1) ^ 2

  calc.equilibrium <- function() {
    # For now works only for g \in R, g >= 2
    coef <- c(-a / b, 1)
    coef[g + 2] <- 1
    coef[is.na(coef)] <- 0
    xs <- solve(polynomial(coef))
    Re(xs[abs(Im(xs)) < 1e-9 & Re(xs) > 0][1])
  }
  calc.eigen <- function() {
    upsilon = -df.star
    values <- df.star * exp(2i * pi / n * (0:(n - 1))) - 1
    circle.points <- df.star * exp(2i * pi / 4 * (0:3)) - 1
    plot <- function() {
      x <- Re(values)
      y <- Im(values)
      labs <- vector("expression", n)
      for (i in 1:n)
        labs[i] <- substitute(expression(lambda[i]), list(i = i - 1))[2]
      graphics::plot(Re(circle.points), Im(circle.points), asp = 1, type = "n", xlab = "", ylab = "", las = 1, main = "Eigenvalues")
      mtext("Re", side = 1, line = 3, las = 1, cex = 1.5)
      mtext("Im", side = 2, line = 2, las = 1, cex = 1.5)
      abline(h = 0)
      abline(v = 0)
      draw.circle(-1, 0, abs(df.star), border = "red", nv = 200)
      points(x, y, pch = 16, col = "blue")
      textxy(x, y, labs, cex = 1.5)
      points(-1, 0, pch = 16, col = "forestgreen")
    }
    list(
      values = values,
      plot = plot
    )
  }
  calc.period <- function(x)
  {
    n <- length(x)
    spec <- spec.ar(c(x), plot = FALSE)
    if (max(spec$spec) > 10) # Arbitrary threshold chosen by trial and error.
    {
      period <- round(1/spec$freq[which.max(spec$spec)])
      if (period == Inf) # Find next local maximum
      {
        j <- which(diff(spec$spec) > 0)
        if (length(j) > 0)
        {
          nextmax <- j[1] + which.max(spec$spec[j[1]:500])
          period <- round(1/spec$freq[nextmax])
        }
        else
          period <- 1
      }
    }
    else
      period <- 1
    return(period)
  }
  calc.repr <- function() {
    x0 <- calc.equilibrium()
    f1 <- df(x0)
    g1 <- -b
    w.star <- sqrt(f1 ^ 2 - b ^ 2)
    T.star <- 2 * pi / w.star
    tau.star <- acos(-g1 / f1) / w.star
    
    simulate <- function(x0 = 1, tau = 1, times = seq(0, 300, by = 0.1), col="red") {
      model <- function(t, x, parms) {
        xp <- ifelse(t < tau, x0, lagvalue(t - tau))
        dx <- f(xp) - x
        list(c(dx), dx = dx, xp = xp)
      }
      names(x0) <- c("x")
      data <- dede(y = x0, times = times, func = model, parms = NULL)
      period <- calc.period(data[-(1:(nrow(data) / 2)), 2]) * (times[2] - times[1])
      plot <- function() {
        df <- data.frame(data)
        ggplot(df, aes(x = time, y = x)) + geom_line(colour = col)
      }
      list(
        x0 = x0,
        tau = tau,
        period = period,
        data = data,
        plot = plot
      )
    }
    
    list(
      x0 = x0,
      f1 = f1,
      g1 = g1,
      w.star = w.star,
      T.star = T.star,
      tau.star = tau.star,
      simulate = simulate
    )
  }

  simulate <- function(x0, times = seq(0, 300, by = 0.1)) {
    model <- function(t, x, params) {
      xp <- c(x[-1], x[1])
      dx <- f(xp) - b * x
      list(dx)
    }
    names(x0) <- paste0("x", 1:length(x0))
    data <- ode(y = x0, times = times, func = model)
    period <- calc.period(data[-(1:(nrow(data) / 2)), 2]) * (times[2] - times[1])
    
    plot.var <- function(col="red") {
      df <- data.frame(data) %>% gather(var, value, -time)
      ggplot(df, aes(x = time, y = value)) + geom_line(colour = col) + facet_grid(. ~ var)
    }
    plot.2d <- function(varx = "x1", vary = "x2", col="red") {
      df <- data.frame(data)
      ggplot(df, aes_string(x = varx, y = vary)) + geom_path(colour = col)
    }
    plot.3d <- function(varx = "x1", vary = "x2", varz = "x3", col="red") {
      lines3d(data[,c(varx, vary, varz)], col = col)
    }
    list(
      data = data, 
      period = period,
      plot.var = plot.var,
      plot.2d = plot.2d,
      plot.3d = plot.3d
    )
  }
  
  print <- function() {
    cat(paste0("f(x) = ", a, "/(1+x^", g, ") - ", ifelse(b == 1, "", b), "x\n"))
    cat(paste0("x0 = ", x.star[1], "\n"))
    cat(paste0("df/dx|x0 = ", df.star[1], "\n"))
    cat(paste0("Repressilator/Hopf: (w = ", repr$w.star, ", T = ", repr$T.star, ", tau = ", repr$tau.star, ")\n"))
    cat(paste0("Eigenvalues = (", paste0(eigen$values, collapse = ", "), ")\n"))
  }
  
  x.star <- rep(calc.equilibrium(), n)
  df.star <- df(x.star)
  repr <- calc.repr()
  eigen <- calc.eigen()
  
  list(
    n = n,
    a = a, 
    g = g, 
    b = b, 
    
    x.star = x.star,
    df.star = df.star,
    repr = repr,
    eigen = eigen,
    
    f = f,
    df = df,
    simulate = simulate,
    print = print
  )
}