library(mlbench)
data("BostonHousing2")

plot(BostonHousing2$town, BostonHousing2$nox)

plot(BostonHousing2$town, BostonHousing2$b)
plot(BostonHousing2$b[BostonHousing2$town=="Boston Downtown"])
