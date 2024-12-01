library(gstat)
library(terra)
library(sf)
library(dplyr)
library(sp)
library(tidyr)
library(parallel)


crs <- "5070"

#coords <-data.frame(expand.grid(x = seq(1, 10000, 10), y =seq(1, 10000, 10)))

x.start <- seq(1, 10000, 10)
y.start <- seq(1, 10000, 10)

field_ind_x <- sort(rep(seq(0,9, 1), 100))
field_ind_y <- sort(rep(seq(0,9, 1), 100))

coords <- data.frame(expand.grid(x = x.start + field_ind_x * 1000, 
                                 y = y.start + field_ind_y * 1000))
coords.test <- filter(coords, x < 6000, y< 6000)



sp::gridded(coords) <- ~x+y
sp::gridded(coords.test) <- ~x+y                     


vgm.model <- vgm(psill = 8, range = 5000, nugget = 1, model = 'Exp')

sim.variable <- function(coords, vgm.model, beta =4){
  
  vg.dummy <- gstat(formula= z~1, locations = ~x+y,  dummy = T,
                    model = vgm.model, nmax = 10, beta = beta)
  
  simmed <- predict(vg.dummy, newdata = coords, nsim = 1)
  
  
  #colnames(data) <- c(c('x', 'y'), paste('var', 1:10, sep = ''))
  
  #write.csv(data, 'sim_geospatial.csv')
  return(simmed)
}

cl <- makeCluster(detectCores())
clusterExport(cl, c('sim.variable', 'coords', 'coords.test'),
              envir = environment())

clusterEvalQ(cl, library(gstat))

fun <- function(vg.mod){sim.variable(coords, vg.mod)$sim1}

vg_mods <- list(
  vgm(psill = 8, range = 100, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 500, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 1000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 1500, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 2000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 3000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 4000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 5000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 6000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 7000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 8000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 9000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 10000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 12000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 14000, nugget = 1, model = 'Exp'),
  vgm(psill = 8, range = 16000, nugget = 1, model = 'Exp')
  
)


fun.test <- function(vg.mod){sim.variable(coords.test, vg.mod)$sim1}

#simmed.vars.test <-  parLapply(cl, vg_mods, fun.test) %>% bind_cols()
#colnames(simmed.vars.test) <- paste('var', 1:16, sep = '')

simmed.vars <- parLapply(cl, vg_mods, fun) %>% bind_cols()
colnames(simmed.vars) <- paste('var', 1:16, sep = '')
data <- simmed.vars
data$x <- coords$x
data$y <- coords$y
data$id <- 1:nrow(data)
write.csv(data, 'sim_geospatial.csv')
