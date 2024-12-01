library(tidyr)
library(dplyr)
library(randomForest)
library(quantregForest)
library(gstat)
library(parallel)

#to do: 
# 1 switch to quantile regression forest and check whether z-scores are z-distributed
#2 investigate variograms of z-transformed prediction errors across x-val strategies.


data <- read.csv('sim_geospatial.csv')

data <- data %>% 
  mutate(target = var1 + var2 + var3 + var4 + exp(var5/10) + log(var6 + abs(min(var6))+ 
              .0001) + var1 * var7 + var3 * exp(var7/10) + var8 * var3 + 
           var9*var5 + var9 * var2 + var10 * var7 + exp(var11/20)*var12 + 
           var13 * var14 * var15 + var16 + var7 * var9 + var1 * exp(var2/20) * var16 +
         var8 * var11)

fl.x <- floor(data$x/1000) + 1
fl.y <- (floor(data$y/1000)) +1
fiel <- fl.x + fl.y *max(fl.x)

field.key <- data.frame(field.1 = unique(fiel), field = rank(unique(fiel)))
data$field.1 <- fiel
data <- merge(data, field.key)
data <- data %>% select(-c(field.1))
data$id <- 1:nrow(data)

always_drop <- c('id', 'x', 'y', 'field', 'target')

sample.fields <- function(data, n.fields, samps.per.field){
  field_ind <- sample(unique(data$field), n.fields)
  
  sample.points <- data.frame()
  for (ind in field_ind){
    sample.points <- rbind(sample.points,
                           sample_n(dplyr::filter(data, field == ind), samps.per.field)
                           )
    
  }
  return (sample.points)
  }


lofo.test <- function(sample_data, masked_vars, paralell = T){
  fun <- function(i){
    test <- dplyr::filter(sample_data, field == i)
    train <- dplyr::filter(sample_data, field != i)
    
    rf <- quantregForest(x = select(train, 
                                    -c(masked_vars, always_drop)), 
                         y =      train$target)
    
    preds_m <- predict(rf, select(test, 
                                  -c(masked_vars, always_drop)),
                       what = mean)
    
    preds_sd <-  predict(rf, select(test, 
                                    -c(masked_vars, always_drop)),
                         what = sd)
    return(data.frame(pred_mean = preds_m, 
                      pred_sd = preds_sd,
                      actual = test$target,
                      x = test$x,
                      y = test$y
                      
    ))  
  }
  
  
  if (paralell){
      cl <- makeCluster(detectCores())
    clusterEvalQ(cl, library(quantregForest))
    clusterEvalQ(cl, library(dplyr))
    clusterExport(cl, c('sample_data', 'masked_vars', 'always_drop'), envir = environment())
    results <- parSapply(cl, unique(sample_data$field), 
                         fun, simplify = F)  %>% bind_rows()
    
  }else{
  results <- data.frame()
  for (ind in unique(sample_data$field)){
   # train <- filter(sample_data, ind != field)
  #  test <- filter(sample_data, ind == field)
    
    #rf <- randomForest( x = select(train, -c(masked_vars, 'x', 'y', 'field', 'target')),
    #                    y = train$target, 
    #                    xtest = select(test, -c(masked_vars, 'x', 'y', 'field', 'target')),
    #                    ytest = test$target
    #                    )
    
    #preds <- rf$test$predicted
   # rf <- quantregForest(x = select(train, 
  #                              -c(masked_vars, 'x', 'y', 'field', 'target')), 
   #                       y =      train$target
    #                      )
    #preds_m <- predict(rf, select(test, 
     #                           -c(masked_vars, 'x', 'y', 'field', 'target')),
      #               what = mean)
    #preds_sd <-  predict(rf, select(test, 
    #                                  -c(masked_vars, 'x', 'y', 'field', 'target')),
    #                       what = sd)
    
    #results <- rbind(results, data.frame(pred_mean = preds_m, 
    #                                     pred_sd = preds_sd,
     #                                    actual = test$target,
     #                                    x = test$x,
    #                                     y = test$y))
    results <- rbind(results, fun(i))
  }
  }
  return(results)
}

losfo.test <- function(sample_data, masked_vars, paralell = T){
  n.samps.per.field <- max(table(sample_data$field))
  sample_data$samp_num <- 1:n.samps.per.field
  
  fun <- function(i){
    test <- dplyr::filter(sample_data, samp_num == i)
    train <- dplyr::filter(sample_data, samp_num != i)
    rf <- quantregForest(x = select(train, 
                                    -c(masked_vars,always_drop, 'samp_num')), 
                         y =      train$target)
    
    preds_m <- predict(rf, select(test, 
                                  -c(masked_vars, always_drop, 'samp_num')),
                       what = mean)
    
    preds_sd <-  predict(rf, select(test, 
                                    -c(masked_vars, always_drop, 'samp_num')),
                         what = sd)
    
    return(data.frame(pred_mean = preds_m, 
                      pred_sd = preds_sd,
                      actual = test$target,
                      x = test$x,
                      y = test$y
                      
    ))  
  }
  
  
  if (paralell){
      cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(quantregForest))
  clusterEvalQ(cl, library(dplyr))
  clusterExport(cl, c('sample_data', 'masked_vars', 'always_drop'), envir = environment())
  results <- parSapply(cl, 1:n.samps.per.field, fun, simplify = F)  %>% bind_rows()
  
  }else{
  results <- data.frame()
  for (i in 1:n.samps.per.field){
    #test <- filter(sample_data, samp_num == i)
    #train <- filter(sample_data, samp_num != i)
    #rf <- randomForest( x = select(train, -c(masked_vars, 'x', 'y', 'field', 'target')),
    #                    y = train$target, 
    #                    xtest = select(test, -c(masked_vars, 'x', 'y', 'field', 'target')),
    #                    ytest = test$target
    #)
    #preds <- rf$test$predicted
    
    
    #rf <- quantregForest(x = select(train, 
    #                                -c(masked_vars, 'x', 'y', 'field', 'target', 'samp_num')), 
    #                                y =      train$target)
    
    #preds_m <- predict(rf, select(test, 
    #                              -c(masked_vars, 'x', 'y', 'field', 'target', 'samp_num')),
    #                   what = mean)
    
    #preds_sd <-  predict(rf, select(test, 
    #                                -c(masked_vars, 'x', 'y', 'field', 'target', 'samp_num')),
    #                     what = sd)
    
    #results <- rbind(results, data.frame(pred_mean = preds_m, 
    #                                     pred_sd = preds_sd,
     #                                    actual = test$target,
    #                                     x = test$x,
    #                                     y = test$y))
    results <- rbind(results, fun(i))
    
  }}
  return(results)
}

rmse <- function(df){
  return(sqrt(mean((df$pred_mean-df$actual )**2)))
}

#'
z_score_var <- function(df){
  return (
    mean(((df$pred_mean - df$actual) / df$pred_sd)**2)
  )
}

get_spe_vg_range <- function(df){
  df$spe <- (df$pred_mean - df$actual) / df$pred_sd
  v <- variogram(spe ~ 1, locations = ~x+y, data = df)
  range <- fit.variogram(v, vgm('Exp'))[2, 'range']
  return(range)
}

compare.tt_splits <- function(full_data, sample_data, masked_vars){
  res_lofo <- lofo.test(sample_data, masked_vars)
  res_losfo <- losfo.test(sample_data, masked_vars)
  
  full_test <- full_data %>% filter(!(id %in% sample_data$id))
  
  full_model <- quantregForest(x = select(sample_data, 
                                  -c(masked_vars, always_drop)), 
                                  y =      sample_data$target
  )
  
  preds_m <- predict(full_model, select(full_test, 
                                -c(masked_vars, 'x', 'y', 'field', 'target')),
                     what = mean)
  
  preds_sd <-  predict(full_model, select(full_test, 
                                  -c(masked_vars, 'x', 'y', 'field', 'target')),
                       what = sd)
  
  out <-  data.frame(pred_mean = preds_m, 
                     pred_sd = preds_sd,
                                       actual = full_test$target, 
                     x = full_test$x, 
                     y = full_test$y)
  
  
  #full_model <- randomForest( x = select(sample_data, -c(masked_vars, 'x', 'y', 'field', 'target')),
  #                                  y = sample_data$target, 
  #                                  xtest = select(full_test, -c(masked_vars, 'x', 'y', 'field', 'target')),
  #                                  ytest = full_test$target)
  
  #preds <- full_model$test$predicted
  
  data_weighted_avg <- left_join(res_lofo, res_losfo%>% select(-actual), 
                             by = c('x', 'y'))
  
  frac_sampled <- length(unique(sample_data$field)) /length(unique(full_data$field))
  
  data_weighted_avg <- data_weighted_avg %>% 
    mutate(pred_mean = pred_mean.x * (1-frac_sampled) + pred_mean.y* frac_sampled,
           pred_sd = sqrt(pred_sd.x**2 * (1-frac_sampled)+ pred_sd.y**2 * frac_sampled))
  
  res <- data.frame(
  
  rmse_lofo = rmse(res_lofo),
  rmse_losfo = rmse(res_losfo),
  rmse_full_set = rmse(out),
  rmse_avged = rmse(data_weighted_avg),
  z_var_lofo = z_score_var(res_lofo),
  z_var_losfo = z_score_var(res_losfo),
  z_var_full = z_score_var(out),
  z_var_avged = z_score_var(data_weighted_avg),
  range_lofo = get_spe_vg_range(res_lofo),
  range_losfo = get_spe_vg_range(res_losfo),
  range_full_set = get_spe_vg_range(sample_n(out, 50000)),
  range_avg = get_spe_vg_range(data_weighted_avg))
  return(res)
  
}

sample_and_test <- function(data, n.fields, samps.per.field){
  sample_data <- sample.fields(data, n.fields, samps.per.field)
  masked_vars <- sample(colnames(data)[1:16] ,5)
  
  compare.tt_splits(data, sample_data, masked_vars)
}


run.experiment <- function(data, n.fields, 
                           samps.per.field, 
                           n.runs 
                           ){
  
  res <- do.call('rbind', sample_and_test(data, n.fields,
                                          samps.per.field
                                                    
                             )) 
  
  return(res%>% mutate(n_fields = n.fields,
                       samps_per_field = samps.per.field))
  }

test.it <- function(){run.experiment(data, 5, 5,10)}

n.fields <- c( 4,5, 10, 20, 25, 50)
samps <- c(250, 125, 100, 50, 25, 20, 10)


res <- data.frame()
for (i in 1:length(n.fields)){
  res <- bind_rows(res, run.experiment(data, n.fields[i], samps[i],
                                   500))
}
write.csv(res, 'results_sim.csv ')

 

res %>% mutate(miss_factor_lofo = rmse_lofo**2 / rmse_full_set**2, 
               miss_factor_losfo = rmse_losfo**2 / rmse_full_set**2,
               miss_factor_mean = ((
  rmse_lofo**2 * (1-(n_fields/100))) + (rmse_losfo**2 * n_fields/100 )) / rmse_full_set**2) %>%
  select(-c(rmse_lofo, rmse_losfo, rmse_full_set)) %>%
  pivot_longer(cols = c(miss_factor_lofo, miss_factor_losfo, miss_factor_mean)) %>%
  #filter(n_fields >5) %>%
  ggplot(aes(y = log(value), color = name)) + geom_boxplot() + 
  geom_hline(aes(yintercept = 0)) + facet_wrap(~n_fields)
                                                                    

res %>% mutate(miss_factor_lofo = rmse_lofo**2 / rmse_full_set**2, 
               miss_factor_losfo = rmse_losfo**2 / rmse_full_set**2,
               miss_factor_mean = ((
                 rmse_lofo**2 * (1-(n_fields/100))) + (rmse_losfo**2 * n_fields/100 )) / rmse_full_set**2) %>%
  select(-c(rmse_lofo, rmse_losfo, rmse_full_set)) %>%
  pivot_longer(cols = c(miss_factor_lofo, miss_factor_losfo, miss_factor_mean)) %>%
  group_by(n_fields, name) %>% summarize(mean(value))

