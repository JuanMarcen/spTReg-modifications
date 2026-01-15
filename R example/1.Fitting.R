# Example of a quantile regression model using convolutions in the coastal covariance matrix

tau <- 0.5

library(qs)
# raw data
df <- qread('R example/data/df_jun_ag.qs')
# stations data
stations <- readRDS('R example/data/stations.rds')

# data preparation
library(dplyr)
library(lubridate)
df <- df %>%
  select(Date, station, Y, l, t, s.1, c.1, elev, dist, g300, g500, g700) %>%
  mutate(
    month = month(Date),
    `t:month6` = ifelse(month == 6, t, 0),
    `t:month7` = ifelse(month == 7, t, 0),
    `t:month8` = ifelse(month == 8, t, 0)
  )


# definition of arguments for MCMC function
Y <- df$Y
X <- matrix(0, nrow = length(Y), ncol = 0) # no fixed effects
V <- cbind(1, df$s.1, df$c.1, df$g300, df$g500, 
           df$g700, df$`t:month6`, df$`t:month7`, 
           df$`t:month8`)

# X_alpha (mean of the GP)
aux.vars <- c('elev', 'dist')
X_alpha <- list()
for (i in 1:ncol(V)){
  X_alpha[[i]] <- cbind(1, unique(df$elev), unique(df$dist))
}

# distance between stations
dist <- readRDS('R example/data/dist.matrix.rds')

# priors of parameters and hyperparameters
M <- matrix(0, nrow = 0, ncol = 1)
P <- matrix(0, nrow = 0, ncol = 0)

M_beta_alpha <- list()
P_beta_alpha <- list()
for (i in 1:length(X_alpha)){
  M_beta_alpha[[i]] <- rep(0, 3)
  P_beta_alpha[[i]] <- 0.0001 * diag(3)
}

da <- 38
db <- 7400
ga <- 2
gb <- 1
ra <- 83
rb <- 24600
na <- 0.1
nb <- 0.1

# initial valkues of parameters and hyperparameters
beta <- matrix(0, nrow = 0, ncol = 1)
alpha <- matrix(0.1, nrow = nrow(stations), ncol = ncol(V))
prec <- 1

hp <- matrix(0, nrow = 5, ncol = ncol(V))
hp[1, ] <- 1 #precision
hp[2, ] <- 3/600 #decay
hp[3, ] <- 1 #varsigma
hp[4, ] <- 3/900 # varphi
hp[5, ] <- 1

beta_alpha <- list()
for (i in 1:length(X_alpha)){
  beta_alpha[[i]] <- rep(0.1, times = ncol(X_alpha[[i]]))
}

# more constants 
N <- nrow(df)
n <- nrow(stations)
p <- 0
r <- ncol(V)

p_alpha <- unlist(lapply(X_alpha, ncol))
s <- rep(0:39, each = 5888)

# number of burning and sampling (as desired)
# to obtain convergence usually a high number
nSims <- 100000
nThin <- 100
nBurnin <- 100000
nReport <- 1000

# more distances
dist_coast <- readRDS('R example/data/dist.vec.rds')
dist_coast_points <- readRDS('R example/data/dist.coast.points2.rds')

dmatcoast_conv <- readRDS('R example/data/phimat.rds')
drmat_conv <- readRDS('R example/data/dr.rds')
lencoast_conv <- drmat_conv[1, 200]

# --- MCMC ---
# call to C++ code
Rcpp::sourceCpp("C++/mcmc.cpp")

# model specification: 1 (no convolution) or 2 (convolution)
model <- 2

# repeat until no error when calculating matrix inverses
repeat {
  conv <- try(spQuantileRcpp(
    tau = tau,
    Y = Y,
    X = X,
    V = V,
    X_alpha = X_alpha,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha,
    P_beta_alpha = P_beta_alpha,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha,
    prec = prec,
    hp = hp,
    beta_alpha = beta_alpha,
    N = N,
    n = n,
    p = p,
    r = r,
    p_alpha = p_alpha,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s,
    model = model
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(conv, "try-error")) {
    return(conv)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}

# traceplots of GP
names <- c('intercept', 's.1', 'c.1', 'g300', 'g500', 'g700', 
           't_month6', 't_month7', 't_month8')

par(mfrow = c(4, 5))
for (j in 0:8){
  for (i in 1:40){
    plot(conv$process[, i + j * 48], type = 'l', main = paste0(names[j + 1], '(', stations$NAME2[i], ')'))
  }
}

# traceplots of parameters in GP
par(mfrow = c(4,2))
for (i in 0:8){
  plot(conv$process[, 41 + i * 48], type = 'l', main = paste0('mu(intercept)_', names[i + 1]))
  plot(conv$process[, 42 + i * 48], type = 'l', main = paste0('mu(elev)_', names[i + 1]))
  plot(conv$process[, 43 + i * 48], type = 'l', main = paste0('mu(dist)_', names[i + 1]))
  plot(conv$process[, 44 + i * 48], type = 'l', main = paste0('1/sigma_k^2_', names[i + 1])) #1/simga_k^2
  plot(conv$process[, 45 + i * 48], type = 'l', main = paste0('decay_', names[i + 1])) # decay
  plot(conv$process[, 46 + i * 48], type = 'l', main = paste0('varsgima_', names[i + 1])) # varsigma
  plot(conv$process[, 47 + i * 48], type = 'l', main = paste0('varphi_', names[i + 1])) #varphi
  plot(conv$process[, 48 + i * 48], type = 'l', main = paste0('prec.coast_', names[i + 1])) # prec.coast
}

# --- KRIGING AND MAP ---
# pre fitted model chains with many iterations
conv <- readRDS('R example/data/conv.trend.q0.50.100k.rds')

# distances for the kriging 
dmatcoast_conv2 <- readRDS('R example/data/phimat.grid.rds')
newcoords <- readRDS('R example/data/grid_km.rds')
coords <- readRDS('R example/data/coords.stations.rds')

dist_coast_points.grid <- readRDS('R example/data/dist.coast.points.grid.rds')
dist_coast.grid <- readRDS('R example/data/dist.vec.grid.rds')

dist_coast_points.comb <- readRDS('R example/data/r.stations.grid.rds')

# function for kriging all GP
kriging <- function(chain, vars){
  
  out <- list()
  
  for (i in 0:(length(vars) - 1)){
    # substract mean of the GP
    w <- chain$process[, 1:40 + i*48]
    Z <- X_alpha[[i + 1]]
    w.def <- w - t(Z %*% t(chain$process[, 41:43 + i*48]))
    
    out[[vars[i + 1]]] <- krigeBayesRcpp(
      w = w.def,
      hp = chain$process[ , 44:48 + i*48, drop = FALSE],
      coords = coords,
      newcoords = newcoords,
      dr = drmat_conv,
      dcoast = dmatcoast_conv,
      newdcoast = dmatcoast_conv2,
      lencoast = lencoast_conv,
      dvec = dist_coast,
      newdvec = dist_coast.grid,
      dmatc = dist_coast_points,
      newdmatc = dist_coast_points.grid,
      combdmatc = dist_coast_points.comb,
      model = model
    )
  }
  
  return(out)
}

vars <- c('intercept', 's.1', 'c.1', 'g300', 'g500', 'g700', 
          't_month6', 't_month7', 't_month8')

# computation time depends on the length of chains
kr.conv <- kriging(conv, vars)

# MAP PLOTS
library(viridis)
library(ggplot2)
library(sf)
library(sp)
limits <- readRDS('R example/data/limits.rds')
background <- readRDS('R example/data/background.rds')
grid <- readRDS('R example/data/grid.rds')

# stations to st object
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")],
      data = stations,
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    ),
    'sf'
  ),
  2062
)

# function to plot GP kriging
plot.gp <- function(kriging, grid, stations, limits, background, var, tau){
  mu <- colMeans(kriging[[var]], na.rm = T)
  
  grid_coords <- cbind(st_coordinates(grid), mu = mu)
  grid_coords <- na.omit(grid_coords)
  grid.new <- st_sf(mu = round(as.vector(mu), 3), geometry = grid)
  
  ggplot(data = background) + 
    geom_sf(fill = "antiquewhite") + 
    xlab("Longitud (º)") + ylab("Latitud (º)") + ggtitle(bquote( .(var) * .(' (') * tau *.(' = ') *.(0.5) *.(')'))) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6, angle = 90),
          axis.title=element_text(size = 10, face = "bold")) + 
    geom_tile(data = grid_coords, aes(X, Y, fill = mu)) +
    geom_tile(data = grid.new, ggplot2::aes(x = st_coordinates(grid.new)[, 1], y = st_coordinates(grid.new)[, 2], fill = mu)) +
    geom_sf(data = stations, aes(color = mu), 
            shape = 21,        # forma con relleno y borde
            color = "black",   # contorno
            size = 3,
            stroke = 1 )+
    scale_fill_gradient2( low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                          space = "Lab", midpoint = 0, limits = c(-5, 5), name = "Distance (km)") +
    scale_color_gradient2( low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                           space = "Lab", midpoint = 0, limits = c(-5, 5), name = "Distance (km)") +
    
    coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])
}

# maps for all GP
# limits in the scale color may be tocuhed for better images
for (var in vars){
  g <- plot.gp(kr.conv, grid, stations, limits, background, var, 0.50)
  print(g)
}



