# FUNCTIONS FOR COMPUTING METRICS IN BAYESIAN MODELS
#fitting of models parallelized
mod_bay <- function(formula, data, tau, vars, 
                    coords, start_beta, inic_procesos,
                    n.samples, n.burnin, n.thin, n.report){
  
  vars <- gsub('`', '', vars)
  ind <- length(vars) + 3
  
  # Bucle de reintento
  repeat {
    mod_try <- try(
      spTm(formula,
           data = data,
           method = 'q',
           quantile = tau,
           coords = coords,
           v = as.matrix(cbind(1, data[, vars])),
           priors = list(
             "beta" = list(M = rep(0, ind), P = 0.0001 * diag(ind)),
             "sigma" = c(0.1, 0.1),
             "phi" = c(38, 7400),
             "mu" = c(0, 0.0001)),
           starting = list(
             "beta" = start_beta,
             "sigma" = 1,
             "alpha" = inic_procesos,
             "hp" = c("mu" = 0, "sigma" = 1, "phi" = 3 / 600)),
           n.samples = n.samples,
           n.burnin = n.burnin,
           n.thin = n.thin,
           n.report = n.report
      ),
      silent = TRUE
    )
    
    # Si NO dio error → salimos del bucle
    if (!inherits(mod_try, "try-error")) {
      return(mod_try)
    }
    
    message("⚠️  Error en spTm() → reintentando...")
    
    # Opcional: pequeña pausa para no saturar la CPU
    Sys.sleep(0.5)
  }
}


traducir_nombres_coef <- function(nombres_coef) { 
  traducidos <- character(length(nombres_coef))
  
  for (i in seq_along(nombres_coef)) {
    nombre <- nombres_coef[i]
    
    if (grepl("^poly\\((.+), 2\\)1$", nombre)) {
      base <- sub("^poly\\((.+), 2\\)1$", "\\1", nombre)
      traducidos[i] <- base
    } else if (grepl("^poly\\((.+), 2\\)2$", nombre)) {
      base <- sub("^poly\\((.+), 2\\)2$", "\\1", nombre)
      traducidos[i] <- paste0("I(", base, "^2)")
    } else {
      traducidos[i] <- nombre
    }
  }
  
  return(traducidos)
}

betas <- function(vars, chain.df){
  vars <- gsub('`', '', vars)
  ind <- length(vars) + 1 #intercept
  params <- chain.df
  
  tr <- traducir_nombres_coef(colnames(params)[2:ind])
  colnames(params)[2:ind] <- gsub('`', '', tr)
  #intercepto
  int <- mean(params[['(Intercept)']])
  int <- rep(int, length=dim(stations)[1])
  
  #fixed_betas
  fixed_betas <- apply(params[, vars], 2, mean)
  fixed_betas <- matrix(rep(fixed_betas, each = dim(stations)[1]), nrow = dim(stations)[1])
  
  elev <- mean(params[,'elev'])
  elev <- matrix(rep(elev, each = dim(stations)[1]), nrow = dim(stations)[1])
  dist <- mean(params[,'dist'])
  dist <- matrix(rep(dist, each = dim(stations)[1]), nrow = dim(stations)[1])
  
  #random_effects (betai(sj))
  cols <- grep('beta', names(params), value=T)
  mu <- apply(params[, cols], 2, mean)
  random_betas <- matrix(mu, nrow = dim(stations)[1])
  
  #put together
  betas <- cbind(int, elev, dist, fixed_betas, random_betas)
  betas <- as.data.frame(betas, row.names = stations$NAME2)
  colnames(betas) <- c('intercept', 'elev', 'dist',
                       vars, paste0('beta',1:(length(vars) + 1)))
  
  
  return(betas)
}

predictions <- function(vars, betas, df, cuantil){
  vars <- gsub('`', '', vars)
  pred <- numeric(nrow(df))
  for (i in 1:dim(stations)[1]){
    ind <- which(df$station == stations$STAID[i])
    for (j in ind){
      #inicializar en interceptos
      comp_esp <- betas[i, 'beta1'] #intercepto espacial para la estacion i
      
      comp_fija <- betas[i, 'intercept'] #intercepto fijo para la estacion i
      
      for (k in 1:(length(vars))){ #beta 1 es la componente espacial del intercepto
        comp_esp <- comp_esp + betas[i, paste0('beta', k + 1)] * df[j, vars[k]]
        comp_fija <- comp_fija + betas[i, vars[k]] * df[j, vars[k]]
      }
      
    
      
      pred[j] <- comp_esp + comp_fija + betas[i, 'elev'] * elev_sc[i] + betas[i, 'dist'] * dist_sc[i]
      
      
    }
  }
  
  return(pred)
}

check <- function(u, tau) {
  return(u * (tau - (u < 0)))  # Implements the quantile loss function
}

R1_bay <- function(pred, tau, data){
  pred_clean <- na.omit(pred)
  
  #dataframe para global
  df <- matrix(NA, nrow=1, ncol=1)
  df <- as.data.frame(df)
  colnames(df) <- c('R1_bay_global')
  
  #dataframe para locales
  df_local <- matrix(NA, nrow=dim(stations)[1], ncol=1)
  df_local <- as.data.frame(df_local, row.names = stations$NAME2)
  colnames(df_local) <- c('R1_bay_local')
  
  #modelos nulos, son para todas variables igual
  mod_nulo_f <- rq(Y ~ as.factor(station), data = data, tau = tau)
  
  rho_estacion <- rep(NA, dim(stations)[1])
  R1_nulo_est <- rep(NA, dim(stations)[1])
  for (j in 1:length(rho_estacion)){
    ind <- which(pred_clean$station == stations$STAID[j])
    rho_estacion[j] <- sum(check(pred_clean$Y[ind] - pred_clean[ind,paste0('pred_q',format(tau, nsmall = 2))],tau = tau))
    R1_nulo_est[j] <- sum(check(mod_nulo_f$residuals[ind], tau = tau))
  }
  
  df['R1_bay_global'] <- 1 - sum(rho_estacion) / mod_nulo_f$rho
  df_local['R1_bay_local'] <- 1 - rho_estacion / R1_nulo_est
  
  return(list(R1_globales = df, R1_locales = df_local))
}

rho_bay <- function(predicciones, tau){
  
  #global
  #dataframe para global
  df <- matrix(NA, nrow=1, ncol=1)
  df <- as.data.frame(df)
  colnames(df) <- c('rho_bay_global')
  pred <- predicciones[[paste0('pred_q',format(tau, nsmall = 2))]]
  dif<- predicciones$Y - pred
  df['rho_bay_global'] <- sum(dif < 0, na.rm = T) / ( 40 * 64 * 92)
  
  #estaciones
  df_est <- matrix(NA, nrow=dim(stations)[1], ncol=1)
  df_est <- as.data.frame(df_est, row.names = stations$NAME2)
  colnames(df_est) <- c('rho_bay_est')
  for (i in 1:dim(stations)[1]){
    ind <- which(predicciones$station == stations$STAID[i])
    dif <- predicciones$Y[ind]-pred[ind]
    df_est[i,] <- sum(dif < 0, na.rm = T) / (64 * 92)
  }
  
  #dias
  df_dia_list <- list() #lista para dias por estacion
  
  day_month <- unique(format(predicciones$Date, "%d-%m"))
  df_dia <- matrix(NA, nrow = length(day_month), ncol=1)
  df_dia <- as.data.frame(df_dia, row.names = day_month)
  colnames(df_dia) <- c('rho_bay_dia')
  
  for (i in 1:length(day_month)){
    ind <- which(format(predicciones$Date, "%d-%m") == day_month[i])
    dif <- predicciones$Y[ind] - pred[ind]
    df_dia[i,] <- sum(dif < 0, na.rm = T) / (64 * 40)
    
    #por estaciones
    for (j in 1:dim(stations)[1]){
      nombre <- stations$NAME2[j]
      
      # Si la estación aún no está en la lista, inicialízala
      if (!(nombre %in% names(df_dia_list))) {
        df_dia_list[[nombre]] <- data.frame(rho_bay_dia = rep(NA, length(day_month)), row.names = day_month)
      }
      df_temp <- df_dia_list[[nombre]]
      
      ind_2 <- which(predicciones$station == stations$STAID[j])
      ind_2 <- ind_2[which(ind_2 %in% ind)]
      dif <- predicciones$Y[ind_2] - pred[ind_2]
      
      #guardado
      df_temp[i, 1] <- sum(dif < 0, na.rm = T) / 64
      df_dia_list[[nombre]] <- df_temp
      
      
    }
    
    
  }
  
  #años
  df_year_list <- list()#lista para años por estacion
  
  year <- unique(year(predicciones$Date))
  df_year <- matrix(NA, nrow = length(year), ncol = 1)
  df_year <- as.data.frame(df_year, row.names = year)
  colnames(df_year) <- c('rho_bay_year')
  for (i in 1:length(year)){
    ind <- which(year(predicciones$Date) == year[i])
    dif <- predicciones$Y[ind] - pred[ind]
    df_year[i,] <- sum(dif < 0, na.rm = T) / (40 * 92)
    
    for (j in 1:dim(stations)[1]){
      nombre <- stations$NAME2[j]
      
      # Si la estación aún no está en la lista, inicialízala
      if (!(nombre %in% names(df_year_list))) {
        df_year_list[[nombre]] <- data.frame(rho_bay_year = rep(NA, length(year)), row.names = year)
      }
      df_temp <- df_year_list[[nombre]]
      
      ind_2 <- which(predicciones$station == stations$STAID[j])
      ind_2 <- ind_2[which(ind_2%in%ind)]
      dif <- predicciones$Y[ind_2] - pred[ind_2]
      
      #guardado
      df_temp[i, 1] <- sum(dif < 0, na.rm = T) / 92
      df_year_list[[nombre]] <- df_temp
    }
    
  }
  
  return(list(rho_globales=df,
              rho_estaciones=df_est,
              rho_años=df_year,
              rho_dias=df_dia,
              rho_dias_est=df_dia_list,
              rho_años_est=df_year_list))
}

#convergence functions
shrink.f.tr.plot <- function(models.list, type, vars){
  L <- length(models.list)
  
  ind <- length(vars) #number of parameters to study
  mod1 <- models.list[[1]]
  mod2 <- models.list[[2]]
  mod3 <- models.list[[3]]
  
  cont <- 0
  if (type == 'intercept'){
    for (i in 1:40){
      chains <- list(
        mod1$p.params.samples[,"(Intercept)"] + 
          mod1$p.params.samples[,"elev"] * scale(elev_sc)[i] + 
          mod1$p.params.samples[,"dist"] * scale(dist_sc)[i] + 
          mod1$p.params.samples[,paste0("beta1(s",i,")")],
        mod2$p.params.samples[,"(Intercept)"] + 
          mod2$p.params.samples[,"elev"] * scale(elev_sc)[i] + 
          mod2$p.params.samples[,"dist"] * scale(dist_sc)[i] + 
          mod2$p.params.samples[,paste0("beta1(s",i,")")],
        mod3$p.params.samples[,"(Intercept)"] + 
          mod3$p.params.samples[,"elev"] * scale(elev_sc)[i] + 
          mod3$p.params.samples[,"dist"] * scale(dist_sc)[i] + 
          mod3$p.params.samples[,paste0("beta1(s",i,")")]
      )
      conv<-gelman.diag(chains,multivariate = F)
    
      cat(paste0("beta1(s",i,")"),": ", round(conv$psrf[1],3), "\n")
      value <- conv$psrf[1] < 1.1
      if (value == TRUE){
        cont <- cont + 1
      }
  
    }
    
    cat('Number of values < 1.1: ', cont,'\n')
    cat('Total number of values:', 40*1, '\n')
    cat('Proportion of convergence: ', cont / 40, '\n')
    
    for (i in 1:40){
      plot(c(mod1$p.params.samples[,"(Intercept)"] + 
               mod1$p.params.samples[,"elev"] * scale(elev_sc)[i] +
               mod1$p.params.samples[,"dist"] * scale(dist_sc)[i] +
               mod1$p.params.samples[,paste0('beta1(s',i,')')]), type = "l",
           xlab='Iteración',ylab=paste0('beta1 + beta1(s',i,')'))
      lines(c(mod2$p.params.samples[,"(Intercept)"] + 
                mod2$p.params.samples[,"elev"] * scale(elev_sc)[i] +
                mod2$p.params.samples[,"dist"] * scale(dist_sc)[i] +
                mod2$p.params.samples[,paste0('beta1(s',i,')')]), col = "gray")
      lines(c(mod3$p.params.samples[,"(Intercept)"] + 
                mod3$p.params.samples[,"elev"] * scale(elev_sc)[i] +
                mod3$p.params.samples[,"dist"] * scale(dist_sc)[i] +
                mod3$p.params.samples[,paste0('beta1(s',i,')')]), col = "salmon")
      
    }
  }else if (type == 'coef'){
    for (j in 2:(ind + 1)){
      cont <- 0
      for (i in 1:40){
        chains <- list(
          mod1$p.params.samples[,j] + 
            mod1$p.params.samples[,paste0("beta",j,"(s",i,")")],
          mod2$p.params.samples[,j] + 
            mod2$p.params.samples[,paste0("beta",j,"(s",i,")")],
          mod3$p.params.samples[,j] + 
            mod3$p.params.samples[,paste0("beta",j,"(s",i,")")]
        )
        conv<-gelman.diag(chains,multivariate = F)
        cat(paste0("beta",j,"(s",i,")"),": ", round(conv$psrf[1],3), "\n")
        
        value <- conv$psrf[1] < 1.1
        if (value == TRUE){
          cont <- cont + 1
        }
      }
      
      cat('Number of values < 1.1: ', cont,'\n')
      cat('Total number of values:', 40*1, '\n')
      cat('Proportion of convergence: ', cont / 40, '\n')
    }
    
    for (j in 2:(ind + 1)){
      for (i in 1:40){
        plot(c(mod1$p.params.samples[,j] + 
                 mod1$p.params.samples[,paste0('beta',j,'(s',i,')')]), type = "l",
             xlab='Iteración',ylab=paste0('beta',j,' + beta',j,'(s',i,')'))
        lines(c(mod2$p.params.samples[,j] + 
                  mod2$p.params.samples[,paste0('beta',j,'(s',i,')')]), col = "gray")
        lines(c(mod3$p.params.samples[,j] + 
                  mod3$p.params.samples[,paste0('beta',j,'(s',i,')')]), col = "salmon")
        
      }
    }
  }
  
}

#ESS
ESS <- function(final.chain, type, vars){
  param_mod_def<-as.mcmc(final.chain)
  
  ind <- length(vars)
  
  if (type == 'intercept'){
    ess_df<-as.data.frame(matrix(NA,ncol=2,nrow=40))
    cont<-1
    for (i in 1:40){
      chains <- list(
        param_mod_def[,"(Intercept)"] + 
          param_mod_def[,"elev"] * scale(elev_sc)[i] + 
          param_mod_def[,"dist"] * scale(dist_sc)[i] + 
          param_mod_def[,paste0("beta1(s",i,")")]
      )
      ess<-effectiveSize(chains)
      ess_df[cont,1]<-paste0("beta1(s",i,")")
      ess_df[cont,2]<-round(ess,3)
      cont<-cont+1
      
    }
    print(ess_df)
    cat('Number of stations with ESS < 200: ',
        dim(ess_df[which(ess_df[,2]<=200),])[1], '\n')
    cat('Station with least ESS: \n')
    print(ess_df[which.min(ess_df[,2]),])
  }else if (type == 'coef'){
    for (j in 2:(ind + 1)){
      ess_df<-as.data.frame(matrix(NA, ncol = 2, nrow = 40))
      cont<-1
      for (i in 1:40){
        chains <- list(
          param_mod_def[,j] + 
            param_mod_def[,paste0("beta",j,"(s",i,")")]
        )
        ess<-effectiveSize(chains)
        ess_df[cont,1]<-paste0("beta",j,"(s",i,")")
        ess_df[cont,2]<-round(ess,3)
        cont<-cont+1
      }
      print(ess_df)
      cat('Number of stations with ESS < 200: ',
          dim(ess_df[which(ess_df[,2]<=200),])[1], '\n')
      cat('Station with least ESS: \n')
      print(ess_df[which.min(ess_df[,2]),])
    }
    
  }
  
}

uncertainty <- function(final.chain, type, vars){
  ind <- length(vars)
  
  if (type == 'intercept'){
    for (i in 1:40){
      q<-quantile(
        final.chain[,"(Intercept)"] + 
          final.chain[,"elev"] * scale(stations_dist$HGHT)[i] + 
          final.chain[,"dist"] * scale(stations_dist$DIST)[i] + 
          final.chain[,paste0("beta1(s",i,")")]
        ,probs=c(0.025,0.975))
      cat(paste0("beta1(s",i,")"),": ", round(q,3), "\n")
    }
  }else if (type == 'coef'){
    for (j in 2:(ind + 1 )){
      ic_df<-as.data.frame(matrix(NA,ncol=3,nrow=40))
      colnames(ic_df)<-c('beta','2.5%','97.5%')
      cont<-1
      for (i in 1:40){
        q <- quantile(
          final.chain[,j] + 
            final.chain[,paste0("beta",j,"(s",i,")")]
          ,probs=c(0.025,0.975))
        ic_df[cont,1]<-paste0("beta",j,"(s",i,")")
        ic_df[cont,2]<-round(q[1],3)
        ic_df[cont,3]<-round(q[2],3)
        cont<-cont+1
      }
      print(ic_df)
    }
  }
  
}
