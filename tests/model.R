library(paratb)
set.seed(42)
set_num_threads(1)
model <- paratb(tspan = seq(1, 10*365, 90))
run(model)
