### working off of Roy's example
# simulation paameters
set.seed(12334) 
sims            <- 100    # how many spatial points to simulate
xlim <- ylim    <- c(0,10) # extent of spatial domain  
res             <- 100     # number of pixels per row/column
mean_samplesize <- 100     # mean N for each binomial sample
posterior_draws <- 100     # how many draws to take from the posterior 



# load (and if needed, install) libraries
for(l in c('RandomFields', 'raster', 'data.table', 'ggplot2', 'lhs', 'INLA', 'gridExtra', 'magick', 'PrevMap', 'geoR')){
  if(!l %in% installed.packages()[,1])
    install.packages(l)
  library(l, character.only = TRUE)
}

# define a helper function to cleanly put table into a raster
insert_raster <- function (raster, new_vals) {
  values(raster)    <- 1
  idx               <- 1:ncell(raster)
  n                 <- ncol(new_vals)
  raster_new        <- raster::brick(replicate(n, raster, simplify = FALSE))
  names(raster_new) <- colnames(new_vals)
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }
  return(raster_new)
}

###############   
#### SIMULATION

# make spatial domain represented as a raster
r <- raster(nrows = res, ncols = res, xmn = xlim[1],ymn = ylim[1], xmx = xlim[2], ymx = ylim[2])


# Simulate an underlying gaussian field
model<-c(var=0.5, range=5, shape=2)
rf    <- RFsimulate(model, x = r, n = 1)
rf_df <- as.data.frame(rf, xy = TRUE) 

# simulate some point locations (s) using latin hypercube sampling
s <- randomLHS(sims, 2) * max(xlim)

# simulate binomial values for those points, get everything into a data.table
d <- data.table(lon = s[,1], lat = s[,2], gp_true = extract(rf, s)) # extract true rf
d[, N      := rpois(.N, mean_samplesize)]  # get a sample size for each point
d[, p_true := plogis( gp_true)]            # true probability surface for comparison
d[, y_obs  := rbinom(.N, N, p_true)]       # simulate the observed value

# plot some of the data
# simulated points
ggplot(d, aes(lon, lat)) + geom_point(aes(size = N), color = 'red') +
  geom_point(aes(size = y_obs)) + theme_minimal()

# simulated points versus true underlying probability
ggplot(d, aes(p_true, y_obs/N)) + geom_point(aes(size = N)) +
  theme_minimal() + geom_abline(intercept=0, slope=1, color='red', lty='dashed')

# true underlying probability surface
ggplot() + geom_raster(data = rf_df, aes(x,y,fill=plogis(sim)))  +
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) + labs(fill = 'Truth') +
  geom_point(data= d, aes(x=lon, y=lat, size = y_obs/N), shape = 1) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank()) 

ggplot() + geom_raster(data = rf_df, aes(x,y,fill=plogis(sim)))  +
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) + labs(fill = 'Truth') +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank())+coord_equal() 

library("geoR")
coords <- as.matrix(d[, c("lon", "lat")])
vari <- variog(coords = coords, data = (d$y_obs/d$N), max.dist=7)
vari.fit <- variofit(vari, ini.cov.pars = c(0.01, 1),
                         cov.model = "matern",
                        fix.nugget = TRUE, nugget = 0 ,
                        fix.kappa = TRUE, kappa = 0.5)
plot(vari)
lines(vari.fit)

c.mcmc <- control.mcmc.MCML(n.sim = 10000, burnin = 2000,
 thin = 8, h = (1.65)/(nrow(d) ^ (1/6)))
par0 <- c(0, vari.fit$cov.pars)

test<-binomial.logistic.MCML(y_obs~ 1, units.m= ~ N, coords= ~lon + lat, 
                       data=d, par0=par0, 
                       start.cov.pars = vari.fit$cov.pars[2], 
                       control.mcmc=c.mcmc,
                         kappa=0.5,
                       fixed.rel.nugget = 0)

par0 <- coef(test)
start <- c(par0[3])


test2<-binomial.logistic.MCML(y_obs~ 1, units.m= ~ N, coords= ~lon + lat, 
                             data=d, par0=par0, 
                             start.cov.pars = start, 
                             control.mcmc=c.mcmc,
                             kappa=0.5,
                             fixed.rel.nugget = 0)

par0 <- coef(test2)
start <- c(par0[3])


test3<-binomial.logistic.MCML(y_obs~ 1, units.m= ~ N, coords= ~lon + lat, 
                             data=d, par0=par0, 
                             start.cov.pars = start, 
                             control.mcmc=c.mcmc,
                             kappa=0.5,
                             fixed.rel.nugget = 0)
par0 <- coef(test3)
start <- c(par0[3])


test4<-binomial.logistic.MCML(y_obs~ 1, units.m= ~ N, coords= ~lon + lat, 
                              data=d, par0=par0, 
                              start.cov.pars = start, 
                              control.mcmc=c.mcmc,
                              kappa=0.5,
                              fixed.rel.nugget = 0)
par0 <- coef(test4)
start <- c(par0[3])


test5<-binomial.logistic.MCML(y_obs~ 1, units.m= ~ N, coords= ~lon + lat, 
                              data=d, par0=par0, 
                              start.cov.pars = start, 
                              control.mcmc=c.mcmc,
                              kappa=0.5,
                              fixed.rel.nugget = 0)



coords <- coordinates(r)
pred.MCML <- spatial.pred.binomial.MCML(test5, coords,
                                        control.mcmc = c.mcmc,
                                        scale.predictions = "prevalence",
                                        standard.errors = TRUE)

crds<-as.data.frame(pred.MCML$grid)
crds$dat<-(pred.MCML$prevalence$predictions)

plot<-ggplot()
plot+geom_tile(data=crds, aes(x, y, fill=dat))+
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) + labs(fill = 'Truth') +

  coord_equal()
  

crds2<-crds[crds$x>3.94 & crds$x<4.0,]
truth<-rf_df[rf_df$x>3.94 & rf_df$x<4.0,]
plot+geom_line(data=truth, aes(1-y, plogis(sim)), size=2, col="red")+geom_line(data=crds2, aes(1-y, dat))+theme_bw()
