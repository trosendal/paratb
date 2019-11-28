## Load the library
library(paratb)

## Initialize the model.
##
## The default model has 10,000 farms with 200 animals in each farm,
## to change the starting state replace u0. In the real model we would
## use a list of herds in the population and the number of animals in
## each at the beginning of the study period:
model <- paratb(tspan = seq(1, 10*365))

## Add 10 infected cows in a herd to start the outbreak
model <- seed_herd(model, 1, "cow", 10)

## Add 10 million random events to the model. This will just inject random
## trade events between herds; in the real model we would use recorded
## animal movements here to reflect the real population.
model <- random_events(model, 1000000)

## Run the model to see the spread
df <- trajectory(run(model))
