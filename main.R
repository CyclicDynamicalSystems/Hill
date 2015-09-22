source("hill-core.R")

h3 <- hill(n = 3, a = 18, g = 3)
h3$print()

h3$eigen$plot()

sim <- h3$simulate(c(1, 2, 3))
sim$period
sim$plot.var()
sim$plot.2d()
sim$plot.3d()

rsim <- h3$repr$simulate(tau = 2.0865227)
rsim$period
rsim$plot()