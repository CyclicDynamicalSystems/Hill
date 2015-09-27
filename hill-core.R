library(polynom)
library(deSolve)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rgl)
library(plotrix)
library(calibrate)
library(xtable)

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
    values <- eigen(equil$A)$values
    vectors <- t(eigen(equil$A)$vectors)
    rownames(vectors) = paste0("e", 1:n)
    colnames(vectors) = paste0("x", 1:n)
    vectors.proj <- vectors
    vectors.proj[seq(1, n, by = 2),] <- Re(vectors.proj[seq(1, n, by = 2),])
    vectors.proj[seq(2, n, by = 2),] <- Im(vectors.proj[seq(2, n, by = 2),])
    vectors.proj <- Re(vectors.proj)

    plot <- function() {
      circle.points <- equil$df[1] * exp(2i * pi / 4 * (0:3)) * 1.2 - 1
      x <- Re(values)
      y <- Im(values)
      labs <- vector("expression", n)
      for (i in 1:n)
        labs[i] <- substitute(expression(lambda[i]), list(i = i))[2]
      graphics::plot(Re(circle.points), Im(circle.points), asp = 1, type = "n", xlab = "", ylab = "", las = 1, main = "Eigenvalues")
      mtext("Re", side = 1, line = 3, las = 1, cex = 1.5)
      mtext("Im", side = 2, line = 2, las = 1, cex = 1.5)
      abline(h = 0)
      abline(v = 0)
      draw.circle(-1, 0, abs(equil$df), border = "red", nv = 200)
      points(x, y, pch = 16, col = "blue")
      textxy(x, y, labs, cex = 1.5)
      points(-1, 0, pch = 16, col = "forestgreen")
    }
    list(
      values = values,
      vectors = vectors,
      vectors.proj = vectors.proj,
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
    equil <- list()
    equil$w <- sqrt(f1 ^ 2 - b ^ 2)
    equil$T <- 2 * pi / equil$w
    equil$tau <- acos(-g1 / f1) / equil$w
    
    simulate <- function(x0 = 1, tau = 1, times = seq(0, 100, by = 0.1), col="red") {
      model <- function(t, x, parms) {
        xp <- ifelse(t < tau, x0, lagvalue(t - tau))
        dx <- f(xp) - x
        list(c(dx), dx = dx, xp = xp)
      }
      names(x0) <- c("x")
      data <- dede(y = x0, times = times, func = model, parms = NULL)
      period <- calc.period(as.numeric(data[-(1:(nrow(data) / 2)), 2])) * (times[2] - times[1])
      plot <- function() {
        df <- data.frame(data)
        ggplot(df, aes(x = time, y = x)) + 
          geom_line(colour = col) +
          ggtitle("Autorepressilator")
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
      equil = equil,
      simulate = simulate
    )
  }

  simulate <- function(x0, times = seq(0, 100, by = 0.1)) {
    model <- function(t, x, params) {
      xp <- c(x[-1], x[1])
      dx <- f(xp) - b * x
      list(dx)
    }
    names(x0) <- paste0("x", 1:length(x0))
    data <- ode(y = x0, times = times, func = model)
    e <- data[,2:(n + 1)] %*% t(eigen$vectors.proj)
    data <- Re(cbind(data, e))
    period <- calc.period(as.numeric(data[-(1:(nrow(data) / 2)), 2])) * (times[2] - times[1])
    
    plot.var <- function(col="red") {
      df <- data.frame(data) %>% gather(var, value, -time)
      ggplot(df, aes(x = time, y = value)) + geom_line(colour = col) + facet_grid(. ~ var)
    }
    plot.2d <- function(varx = "x1", vary = "x2", col="red") {
      df <- data.frame(data)
      ggplot(df, aes_string(x = varx, y = vary)) + 
        geom_path(colour = col) +
        ggtitle(paste0("Phase portrait (", varx, ", ", vary, ")"))
    }
    plot.3d <- function(varx = "x1", vary = "x2", varz = "x3", col="red") {
      open3d()
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
  
  simulate.multi <- function(size = 3, lim = 5, times = seq(0, 300, by = 0.1)) {
    data <- list()
    for (i in 1:size)
      data[[i]] <- simulate(runif(n) * lim, times)
    plot.2d <- function(varx = "x1", vary = "x2") {
      df <- do.call(rbind, 
        lapply(seq_along(data), function(i) data.frame(data[[i]]$data) %>% mutate(id = as.factor(paste0("t", i)))))
      ggplot(df, aes_string(x = varx, y = vary)) + 
        geom_path(aes(colour = id)) +
        ggtitle(paste0("Phase portrait (", varx, ", ", vary, ")"))
    }
    plot.3d <- function(varx = "x1", vary = "x2", varz = "x3") {
      open3d()
      cols <- rainbow(size)
      for (i in 1:size)
        lines3d(data[[i]]$data[,c(varx, vary, varz)], col = cols[i])
    }
    list(
      data = data,
      size = size,
      plot.2d = plot.2d,
      plot.3d = plot.3d
    )
  }
  
  overview <- function() {
    mat <- function(caption, matrix) {
      paste0(caption, ":\n", paste0(paste0("  ", capture.output(matrix)), collapse = "\n"), "\n\n")
    }
    val <- function(caption, value) {
      paste0("  ", caption, " = ", value, "\n")
    }
    paste0(
      paste0("f(x) = ", a, "/(1+x^", g, ") - ", ifelse(b == 1, "", b), "x\n"),
      paste0("\n"),
      val("x*", equil$x[1]),
      val("df/dx|x*", equil$df[1]),
      paste0("\n"),
      mat("Jacobian matrix", equil$A),
      paste0("Eigenvalues:\n", paste0("  lambda_", 1:n, " = ", eigen$values, collapse = "\n"), "\n\n"),
      mat("Eigenvectors", eigen$vectors),
      paste0("Repressilator:\n"),
      val("w*", repr$equil$w),
      val("T*", repr$equil$T),
      val("tau*", repr$equil$tau),
      paste0("")
    )
  }
  
  print <- function() {
    cat(overview())
  }
  
  equil <- list()
  equil$x <- rep(calc.equilibrium(), n)
  equil$df <- df(equil$x)
  equil$A <- matrix(data = rep(0, n * n), nrow = n, ncol = n)
  diag(equil$A) <- -1
  equil$A[(row(equil$A) - col(equil$A) + n) %% n == 1] <- equil$df
  
  repr <- calc.repr()
  eigen <- calc.eigen()
  
  list(
    n = n,
    a = a, 
    g = g, 
    b = b, 
    
    equil = equil,
    repr = repr,
    eigen = eigen,
    
    f = f,
    df = df,
    simulate = simulate,
    simulate.multi = simulate.multi,
    print = print,
    overview = overview
  )
}