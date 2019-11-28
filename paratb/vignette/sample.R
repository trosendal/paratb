library(paratb)

model<- seed_herd(paratb(tspan = 1:1000), 1, "cow", 10)
model<- seed_herd(model, 2, "cow", 10)
model<- seed_herd(model, 3, "cow", 10)
model<- seed_herd(model, 4, "cow", 10)

paratb::counti(run(model))
df <- trajectory(run(model))

df[df$time == 300,]
