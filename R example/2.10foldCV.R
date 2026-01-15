# 10-fold CROSS-VALIDATION
rm(list = ls())

tau <- 0.5

library(qs)
df <- qread('R example/data/df_jun_ag.qs')
stations <- readRDS('R example/data/stations.rds')
Rcpp::sourceCpp("C++/mcmc.cpp")

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


# fucntions to calcualte R1 (and more)
source('R example/src/metrics_bay.R')
library(quantreg)

set.seed(05052002)
# Folds creation
stations$group <- sample(rep(1:10, each = 4))
stations$r <- readRDS('R example/data/r.stations.rds')

#R1.CV <- data.frame(matrix(NA, ncol = 1, nrow = 10))
dist <- readRDS('R example/data/dist.matrix.rds')
dist_coast <- readRDS('R example/data/dist.vec.rds')
dist_coast_points <- readRDS('R example/data/dist.coast.points2.rds')
dmatcoast_conv <- readRDS('R example/data/phimat.rds')
drmat_conv <- readRDS('R example/data/dr.rds')
dmatcoast_conv2 <- readRDS('R example/data/phimat.grid.rds')
coords <- readRDS('R example/data/coords.stations.rds')

cv.fold <- function(j, df, stations, tau, model){
  test <- which(stations$group == j)
  
  #rownames(R1.CV)[j] <- paste(test, collapse = '-')
  
  train <- setdiff(1:40, test)
  
  test.df <- df[which(df$station %in% stations$STAID[test]), ]
  train.df <- df[-which(df$station %in% stations$STAID[test]), ] 
  
  stations.train <- stations[train, ]
  stations.test <- stations[test, ]
  
  # definition of arguments for MCMC function
  Y <- train.df$Y
  X <- matrix(0, nrow = length(Y), ncol = 0)
  V <- cbind(1, train.df$s.1, train.df$c.1, train.df$g300, train.df$g500, train.df$g700, train.df$`t:month6`, train.df$`t:month7`, train.df$`t:month8`)
  
  # X_alpha
  aux.vars <- c('elev', 'dist')
  X_alpha <- list()
  for (i in 1:ncol(V)){
    X_alpha[[i]] <- cbind(1, unique(train.df$elev), unique(train.df$dist))
  }
  
  X_alpha.test <- list()
  for (i in 1:ncol(V)){
    X_alpha.test[[i]] <- cbind(1, unique(test.df$elev), unique(test.df$dist))
  }
  
  
  dist.train <- dist[train, train]
  dist.test <- dist[test, test]
  
  #priors
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
  
  # initial
  beta <- matrix(0, nrow = 0, ncol = 1)
  alpha <- matrix(0.1, nrow = nrow(stations.train), ncol = ncol(V))
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
  N <- nrow(train.df)
  n <- nrow(stations.train)
  p <- 0
  r <- ncol(V)
  
  p_alpha <- unlist(lapply(X_alpha, ncol))
  s <- rep(0:35 , each = 5888)
  
  nSims <- 10
  nThin <- 1
  nBurnin <- 10
  nReport <- 1
  
  #more distances
  dist_coast.train <- dist_coast[train]
  dist_coast.test <- dist_coast[test]
  dist_coast_points.train <- dist_coast_points[train, train]
  dist_coast_points.test <- dist_coast_points[test, test]
  
  
  dmatcoast_conv.train <- dmatcoast_conv[train, ]
  dmatcoast_conv.test <- dmatcoast_conv[test, ]
  
  
  lencoast_conv <- drmat_conv[1, 200]
  
  # TRAINING
  repeat {
    basura <- try(spQuantileRcpp(
      tau = tau,
      Y = Y,
      X = X,
      V = V,
      X_alpha = X_alpha,
      dist = dist.train,
      dist_coast = dist_coast.train,
      dist_coast_point = dist_coast_points.train,
      dmatcoast_conv = dmatcoast_conv.train,
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
    if (!inherits(basura, "try-error")) {
      break
    }
    
    message("⚠️  Error en spTm() → reintentando...")
    
    # Opcional: pequeña pausa para no saturar la CPU
    Sys.sleep(0.5)
  }
  
  # prediction in new stations
  # for launching in big computer
  
  # distances
  
  coords.train <- coords[train, ]
  coords.test <- coords[test, ]
  
  
  
  # cambiar para cuando no sea conv !!!! 36 x 4
  dist_coast_points.comb <- abs(outer(stations.train$r, stations.test$r, '-'))
  
  # function for kriging all parameters
  kriging <- function(chain, vars){
    
    out <- list()
    
    for (i in 0:(length(vars) - 1)){
      # substract mean of the GP
      w <- chain$process[, 1:36 + i*44]
      Z <- X_alpha[[i + 1]]
      w.def <- w - t(Z %*% t(chain$process[, 37:39 + i*44]))
      
      out[[vars[i + 1]]] <- krigeBayesRcpp(
        w = w.def,
        hp = chain$process[ , 40:44 + i*44, drop = FALSE],
        coords = coords.train,
        newcoords = coords.test,
        dr = drmat_conv,
        dcoast = dmatcoast_conv.train,
        newdcoast = dmatcoast_conv.test,
        lencoast = lencoast_conv,
        dvec = dist_coast.train,
        newdvec = dist_coast.test,
        dmatc = dist_coast_points.train,
        newdmatc = dist_coast_points.test,
        combdmatc = dist_coast_points.comb,
        model = model
      )
      
      #add mean of GP for predictions
      Z <- X_alpha.test[[i + 1]]
      out[[vars[i + 1]]] <- out[[vars[i + 1]]] + t(Z %*% t(basura$process[, 37:39 + i*44]))
      
    }
    
    return(out)
  }
  
  vars <- c('intercept', 's.1', 'c.1', 'g300', 'g500', 'g700', 
            't_month6', 't_month7', 't_month8')
  
  basura2 <- kriging(basura, vars)
  
  betas <- do.call(cbind, basura2)
  betas <- matrix(colMeans(betas))
  
  pred <- numeric(nrow(test.df))
  for (i in 1:4){
    ind <- which(test.df$station == stations.test$STAID[i])
    X <- as.matrix(cbind(1, test.df[ind, c('s.1', 'c.1', 'g300', 'g500', 'g700', 
                                           't:month6', 't:month7', 't:month8')]))
    betas.aux <- betas[c(i, i + 4, i + 4 * 2, i + 4 * 3, i + 4 * 4, i + 4 * 5,
                         i + 4 * 6, i + 4 * 7, i + 4 * 8), ,drop = FALSE]
    pred[ind] <- X %*% betas.aux
  }
  
  pred <- cbind(test.df[,c('Date', 'station', 'Y')], pred)
  colnames(pred)[4] <- 'pred_q0.50'
  
  
  R1 <- R1_bay(pred, 0.50, test.df)
  
  out <- mean(R1$R1_locales[, 1], na.rm = T) # what I save
  return(list(
    fold = j,
    test_ids = paste(test, collapse = "-"),
    R1.mean = out,
    R1 = na.omit(R1$R1_locales),
    model = basura,
    kriging = basura2
  )
  )
}


# COMPUTATION IN PARALLEL
library(parallel)

ncores <- min(10, detectCores() - 1)
cl <- makeCluster(ncores)

# Export variables to the clusters
clusterExport(
  cl,
  varlist = c(
    "df", "stations", "tau", "dist", "dist_coast", "dist_coast_points",
    "dmatcoast_conv", "dmatcoast_conv2", "drmat_conv","coords",
    "R1_bay", "check"
  ),
  envir = environment()
)

# Export c++ code and needed libraries to clusters
clusterEvalQ(cl, {
  library(Rcpp)
  Rcpp::sourceCpp("C++/mcmc.cpp")
  library(quantreg)
})

# model with no convolution
res.coastal <- parLapply(cl, 1:10, cv.fold,
                         df = df,
                         stations = stations,
                         tau = tau, 
                         model = 1)

# model with convolution
res.conv <- parLapply(cl, 1:10, cv.fold,
                      df = df,
                      stations = stations,
                      tau = tau, 
                      model = 2)

# stop the parallelization
stopCluster(cl)

# results that may be compared
R1.coastal <- data.frame(
  R1 = sapply(res.coastal, `[[`, "R1.mean"),
  row.names = sapply(res.coastal, `[[`, "test_ids")
)
R1.conv <- data.frame(
  R1 = sapply(res.conv, `[[`, "R1.mean"),
  row.names = sapply(res.conv, `[[`, "test_ids")
)

res <- data.frame(conv = R1.conv, coastal = R1.coastal)
names(res) <- c('R1.conv', 'R1.coastal')
apply(res, 2, mean)
