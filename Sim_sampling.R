library(tidyr)
library(dplyr)
library(randomForest)
library(quantregForest)
library(gstat)
library(paralell)

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
data <- data %>% select(-c(fiel))


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


lofo.test <- function(sample_data, masked_vars){
  
  results <- data.frame()
  for (ind in unique(sample_data$field)){
    train <- filter(sample_data, ind != field)
    test <- filter(sample_data, ind == field)
    
    rf <- randomForest( x = select(train, -c(masked_vars, 'x', 'y', 'field', 'target')),
                        y = train$target, 
                        xtest = select(test, -c(masked_vars, 'x', 'y', 'field', 'target')),
                        ytest = test$target
                        )
    preds <- rf$test$predicted
    
    results <- rbind(results, data.frame(pred = preds, actual = test$target))
  }
    
  return(results)
}

losfo.test <- function(sample_data, masked_vars){
  n.samps.per.field <- max(table(sample_data$field))
  sample_data$samp_num <- 1:n.samps.per.field
  results <- data.frame()
  for (i in 1:n.samps.per.field){
    test <- filter(sample_data, samp_num == i)
    train <- filter(sample_data, samp_num != i)
    rf <- randomForest( x = select(train, -c(masked_vars, 'x', 'y', 'field', 'target')),
                        y = train$target, 
                        xtest = select(test, -c(masked_vars, 'x', 'y', 'field', 'target')),
                        ytest = test$target
    )
    preds <- rf$test$predicted
    
    results <- rbind(results, data.frame(pred = preds, actual = test$target))
  }
  return(results)
}

rmse <- function(df){
  return(sqrt(mean((df$pred-df$actual )**2)))
}

compare.tt_splits <- function(full_data, sample_data, masked_vars){
  res_lofo <- lofo.test(sample_data, masked_vars)
  res_losfo <- losfo.test(sample_data, masked_vars)
  
  full_test <- full_data %>% filter(!(id %in% sample_data$id))
  full_model <- randomForest( x = select(sample_data, -c(masked_vars, 'x', 'y', 'field', 'target')),
                                    y = sample_data$target, 
                                    xtest = select(full_test, -c(masked_vars, 'x', 'y', 'field', 'target')),
                                    ytest = full_test$target)
  
  preds <- full_model$test$predicted
  
  out <- data.frame(pred = preds, actual = full_test$target)
  
  
  res <- data.frame(
  rmse_lofo = rmse(res_lofo),
  rmse_losfo = rmse(res_losfo),
  rmse_full_set = rmse(out))
  return(res)
  
}

sample_and_test <- function(data, n.fields, samps.per.field){
  sample_data <- sample.fields(data, n.fields, samps.per.field)
  masked_vars <- sample(colnames(data)[1:16] ,5)
  
  compare.tt_splits(data, sample_data, masked_vars)
}


run.experiment <- function(data, n.fields, 
                           samps.per.field, 
                           n.runs, 
                           paralell = T){
  fun <- function(i){sample_and_test(data, n.fields,samps.per.field)}
  if (paralell){
    fun <- function(i){sample_and_test(data, n.fields,
                                       samps.per.field)}
    cl <- makeCluster(detectCores())
    clusterExport(cl, c('data', 'sample_and_test', 'sample.fields', 'compare.tt_splits',
                        'lofo.test', 'losfo.test',
                        'n.fields', 'samps.per.field', 'rmse' ),
                  envir = environment())
    
    clusterEvalQ(cl, library(gstat))
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(randomForest))
    
    res <- parLapply(cl, as.list(1:n.runs), fun) %>% bind_rows()
  }else{
  res <- do.call('rbind', replicate(n.runs, fun(1)
                                                    ,
                             simplify = FALSE)) 
  }
  return(res%>% mutate(n_fields = n.fields,
                       samps_per_field = samps.per.field))
  }

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

