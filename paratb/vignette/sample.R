## Load the library
library(paratb)

## Initialize the model.
##
## The default model has 10,000 farms with 200 animals in each farm,
## to change the starting state replace u0. In the real model we would
## use a list of herds in the population and the number of animals in
## each at the beginning of the study period:
model <- paratb(tspan = seq(1, 10*365, 90))

## Seed infection in 0.2% of the herds
for(i in seq_len(ncol(model@u0)*0.002)) {
model <- seed_herd(model,
                   i,
                   "cow",
                   ceiling(model@u0["S3",i] * 0.02))
model <- seed_herd(model,
                   i,
                   "heifer",
                   ceiling(model@u0["S2",i] * 0.02))
model <- seed_herd(model,
                   i,
                   "calf",
                   ceiling(model@u0["S1",i] * 0.02))
}

## Add 5 million random events to the model. This will just inject random
## trade events between herds; in the real model we would use recorded
## animal movements here to reflect the real population.
model <- random_events(model, 1000000)

## Run the model
result <- run(model)

## Just extract the data we need about the herd status at each time point
df <- counti(result, "true")

plot(df$time,
     df$count,
     type = "s",
     ylim = c(0, max(df$count)),
     ylab = "Count of positive herds",
     xlab = "Time in days")


## Now run the same model but inject bulk milk sampling events in 1000
## herds per year and monitor the outcome of that compared to the
## disease.

years <- max(model@tspan)/365

herds_for_surveillance <- sample(seq_len(ncol(model@u0)),
                                 1000 * years,
                                 replace = TRUE)

times_for_surveillance <- sample(min(model@tspan):max(model@tspan),
                                 length(herds_for_surveillance),
                                 replace = TRUE)

model <- add_surveillance_event(model,
                                i = herds_for_surveillance,
                                t = times_for_surveillance,
                                n = 1L)

## Run the model again but with the surveillance events
result <- run(model)

## Extract the data we need about the herd status detection status
df <- counti(result, "true")
df_bulk_milk <- counti(result, "milk")

plot(df$time,
     df$count,
     type = "s",
     ylim = c(0, max(df$count)),
     ylab = "Count of positive herds",
     xlab = "Time in days")

lines(df_bulk_milk$time,
     df_bulk_milk$count,
     type = "s",
     col = "#939393",
     lty = 2)
