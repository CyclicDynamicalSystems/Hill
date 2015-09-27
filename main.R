source("hill-core.R")

h3 <- hill(n = 3, a = 18, g = 3)
h3$print()

h3$eigen$plot()

sim <- h3$simulate(c(1, 2, 3))
sim$plot.2d("e2", "e3")
sim$plot.3d("e1", "e2", "e3")
sim$period
sim$plot.var()
sim$plot.2d()
sim$plot.3d()

rsim <- h3$repr$simulate(tau = 2.0865227)
rsim$period
rsim$plot()

sim10 <- h3$simulate.multi(10)
sim10$plot.2d("e2", "e3")
sim10$plot.3d("e1", "e2", "e3")
