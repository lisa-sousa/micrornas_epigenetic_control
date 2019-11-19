
##############################
#Libraries                   #
##############################

library(glmnet)
library(Cairo)
library(matrixStats)
library(RColorBrewer)
colors = brewer.pal(9,'Set1')

##############################
#READ DATA                   #
##############################

load('/project/histone_marks_miRNA/data/RData_classification/hela_pre_miRNA_promoter_class.RData')
head(data_set,4)

output_dir = '/project/histone_marks_miRNA/pipeline/classification/plots/'
pdf_sub_name = 'pre_miRNA_promoter_new'

###################################
#REGULARIZED LOGISTIC REGRESSION: #
###################################

####################################
#Prediction with test/training set:
####################################

#parameters
runs = 500
sample_size_training_data = c(1000,1000) #for target: 0/1
sample_size_cv_data = c(800,800) #for target: 0/1
alpha = 0.5 #1:Lasso, 0.5:ElasticNet, 0:Ridge
if(alpha == 0){reg = 'ridge'}; if(alpha == 0.5){reg = 'elnet'}; if(alpha == 1){reg = 'lasso'}
thr_class_dist = 0.25
thr_coef = 0.05

#setup and run regularized logistic regression
counter = 1
coef_variance = as.data.frame(matrix(ncol=ncol(data_set),nrow=runs))
colnames(coef_variance) = colnames(data_set)
true_pos_rate_variance = c()
true_neg_rate_variance = c()

while(counter < runs){
  print(counter)
  best_lambdas = c()
  best_accuracy = c()

  #get training/test set
  target_index_0 = grep(0,target)
  target_index_1 = grep(1,target)
  training_index = sort(c(sample(x=target_index_0, size=sample_size_training_data[1], replace=F),
                          sample(x=target_index_1, size=sample_size_training_data[2], replace=F)))
  
  training_target = target[training_index]
  training_data_set = data_set[training_index,]
  
  test_target = target[-training_index]
  test_data_set = data_set[-training_index,]
  
  #run model 200 times to get sequence of optimal lambdas with manuel cv
  while(length(best_lambdas) < 200){
    training_target_index_0 = grep(0,training_target)
    training_target_index_1 = grep(1,training_target)
    cv_index = sort(c(sample(x=training_target_index_0, size=sample_size_cv_data[1], replace=F),
                      sample(x=training_target_index_1, size=sample_size_cv_data[2], replace=F)))
    
    cv_data = training_data_set[cv_index,]
    cv_label = training_target[cv_index]
    validation_data = training_data_set[-cv_index,]
    validation_label = training_target[-cv_index]
    
    lmodel_binom = glmnet(x = data.matrix(cv_data), y = cv_label, family='binomial',standardize=F,nlambda=100,alpha=alpha) 
    predictions = predict(lmodel_binom, data.matrix(validation_data), type="class")
    
    accuracy = c()
    distance = c()
    lambda = lmodel_binom$lambda
    
    for(p in 1:ncol(predictions)){
      fit <- data.frame(true_label = validation_label,predicted_label = as.factor(predictions[,p]))
      
      if(sum(fit$predicted_label == 1) == 0){
        TP = table(fit)[1,1]; FP = table(fit)[2,1]; TN = 0; FN = 0
      }else if(sum(fit$predicted_label == 0) == 0){
        TP = 0; FP = 0; TN = table(fit)[2,1]; FN = table(fit)[1,1]
      }else{
        TP = table(fit)[1,1]; FP = table(fit)[2,1]; TN = table(fit)[2,2]; FN = table(fit)[1,2]
      }
      
      true_pos_rate = TP / (TP + FN)
      true_neg_rate = TN / (TN + FP)
      accuracy = c(accuracy,mean(c(true_pos_rate,true_neg_rate)))
      distance = c(distance,round(abs(true_pos_rate-true_neg_rate),2))
    }
    reduced_accuracy = accuracy[distance < thr_class_dist]
    reduced_lambda = lambda[distance < thr_class_dist]
    best_lambdas = c(best_lambdas,reduced_lambda[which.max(reduced_accuracy)]) 
    best_accuracy = c(best_accuracy,reduced_accuracy[which.max(reduced_accuracy)])
  }
  
  training_target_index_0 = grep(0,training_target)
  training_target_index_1 = grep(1,training_target)
  cv_index = sort(c(sample(x=training_target_index_0, size=sample_size_cv_data[1], replace=F),
                    sample(x=training_target_index_1, size=sample_size_cv_data[2], replace=F)))
  
  cv_data = training_data_set[cv_index,]
  cv_label = training_target[cv_index]
  validation_data = training_data_set[-cv_index,]
  validation_label = training_target[-cv_index]
  
  lmodel_binom = glmnet(x = data.matrix(cv_data), y = cv_label, family='binomial',standardize=F,lambda=best_lambdas,alpha=alpha) 
  predictions = predict(lmodel_binom, data.matrix(validation_data), type="class")
  
  accuracy = c()
  distance = c()
  lambda = lmodel_binom$lambda
  
  for(p in 1:ncol(predictions)){
    fit <- data.frame(true_label = validation_label,predicted_label = as.factor(predictions[,p]))
    
    if(sum(fit$predicted_label == 1) == 0){
      TP = table(fit)[1,1]; FP = table(fit)[2,1]; TN = 0; FN = 0
    }else if(sum(fit$predicted_label == 0) == 0){
      TP = 0; FP = 0; TN = table(fit)[2,1]; FN = table(fit)[1,1]
    }else{
      TP = table(fit)[1,1]; FP = table(fit)[2,1]; TN = table(fit)[2,2]; FN = table(fit)[1,2]
    }
    
    true_pos_rate = TP / (TP + FN)
    true_neg_rate = TN / (TN + FP)
    accuracy = c(accuracy,mean(c(true_pos_rate,true_neg_rate)))
    distance = c(distance,round(abs(true_pos_rate-true_neg_rate),2))
  }
  reduced_accuracy = accuracy[distance < thr_class_dist]
  reduced_lambda = lambda[distance < thr_class_dist]
  best_lambda = reduced_lambda[which.max(reduced_accuracy)]
  
  
  #predicting on test set
  predictions = predict(lmodel_binom, data.matrix(test_data_set), type="class", s=best_lambda)
  
  accuracy = c()
  distance = c()
  for(p in 1:ncol(predictions)){
    fit <- data.frame(true_label = test_target,predicted_label = as.factor(predictions[,p]))
    
    if(sum(fit$predicted_label == 1) == 0){
      TP = table(fit)[1,1]; FP = table(fit)[2,1]; TN = 0; FN = 0
    }else if(sum(fit$predicted_label == 0) == 0){
      TP = 0; FP = 0; TN = table(fit)[2,1]; FN = table(fit)[1,1]
    }else{
      TP = table(fit)[1,1]; FP = table(fit)[2,1]; TN = table(fit)[2,2]; FN = table(fit)[1,2]
    }
    
    true_pos_rate = TP / (TP + FN)
    true_neg_rate = TN / (TN + FP)
    accuracy = c(accuracy,mean(c(true_pos_rate,true_neg_rate)))
    distance = c(distance,round(abs(true_pos_rate-true_neg_rate),2))
  }
  
  reduced_accuracy = accuracy[distance < thr_class_dist]
  if(length(reduced_accuracy) != 0){
    reduced_predictions = as.data.frame(predictions[,distance < thr_class_dist])
    reduced_coefs = as.data.frame(as.matrix(lmodel_binom$beta[,distance < thr_class_dist]))
    counter = counter + 1
    best_prediction = reduced_predictions[,which.max(reduced_accuracy)]
    coef_variance[counter,] = reduced_coefs[,which.max(reduced_accuracy)]
    best_fit <- data.frame(true_label = test_target,predicted_label = as.factor(best_prediction))
    
    TP = table(best_fit)[1,1]; FP = table(best_fit)[2,1]; TN = table(best_fit)[2,2]; FN = table(best_fit)[1,2]
    true_pos_rate = TP / (TP + FN)
    true_neg_rate = TN / (TN + FP)
    print('accuracy')
    print(mean(c(true_neg_rate,true_pos_rate)))
    
    true_pos_rate_variance = c(true_pos_rate_variance,true_pos_rate)
    true_neg_rate_variance = c(true_neg_rate_variance,true_neg_rate)
  }else{print('Criterium for class accuracy not met')}
}


#plot variance in predictions
accuracy_df = data.frame(true_pos_rate = true_pos_rate_variance,true_neg_rate = true_neg_rate_variance)
accuracy_df = cbind(accuracy_df,rowMeans(accuracy_df))
colnames(accuracy_df) = c('TP_rate','TN_rate','accuracy')

coef_variance = na.omit(coef_variance)
coef_variance_sorted = coef_variance

#shrinking outliers for plot
for(i in 1:ncol(coef_variance_sorted)){
  column = coef_variance_sorted[,i]
  b = boxplot(column)
  if(length(b$out)>0){
    for(j in 1:length(b$out)){
      if(b$out[j] > b$stats[5,1]){column[column==b$out[j]]=b$stats[5,1]}
      if(b$out[j] < b$stats[1,1]){column[column==b$out[j]]=b$stats[1,1]}
    }
  }
  coef_variance_sorted[,i] = column
}

#sort by decreasing median
median_sort = colMedians(as.matrix(coef_variance_sorted))
names(median_sort) = colnames(coef_variance_sorted)
coef_variance_sorted_median = coef_variance_sorted[,c(names(sort(median_sort, decreasing=TRUE)))]

#if table with plotting names exist use them
if(exists('class_table')){
  old_names = names(coef_variance_sorted_median)
  new_names = c()
  for(name in old_names){
    new_name = class_table$abbr[class_table$abbr_GSE == name]
    if(length(new_name) == 0){
      new_name = name
    }
    new_names = append(new_names,new_name)
  }
  names(coef_variance_sorted_median) = new_names
  head(coef_variance_sorted_median,4)
}

#sort by decreasing mean
mean_sort = colMeans(coef_variance_sorted)
names(mean_sort) = colnames(coef_variance_sorted)
coef_variance_sorted_mean = coef_variance_sorted[,c(names(sort(mean_sort, decreasing=TRUE)))]

#if table with plotting names exist use them
if(exists('class_table')){
  old_names = names(coef_variance_sorted_mean)
  new_names = c()
  for(name in old_names){
    new_name = class_table$abbr[class_table$abbr_GSE == name]
    if(length(new_name) == 0){
      new_name = name
    }
    new_names = append(new_names,new_name)
  }
  names(coef_variance_sorted_mean) = new_names
  head(coef_variance_sorted_mean,4)
}

#get mean coefficients and remove features with coefficients smaller than thr_coef
coefs_vec = colMeans(coef_variance_sorted_mean)
names(coefs_vec) = colnames(coef_variance_sorted_mean)
reduced_vector = coefs_vec[abs(coefs_vec) > thr_coef] 

#plot
estWidth <- min(4 * log(ncol(data_set)), 4 * log(nrow(data_set)))
estHeight <- min(4 * log(ncol(data_set)), 4 * log(nrow(data_set)))
CairoPDF(file = paste(output_dir,reg,'_regularized_logistic_regression_training_test_',pdf_sub_name,'.pdf',sep=''), 
         width = estWidth, height = estHeight)

plot(1:nrow(accuracy_df),accuracy_df$accuracy,ylab='accuracy',xlab='run',type='l',col=colors[1],ylim=c(0,1),main='Accuracy')
lines(1:nrow(accuracy_df),accuracy_df$TP_rate,col=colors[2])
lines(1:nrow(accuracy_df),accuracy_df$TN_rate,col=colors[3])
legend ("topright", c("total", "class1", "class2",paste('regularization:',reg)), col=c(colors[1:3],'black'), lty=1, lwd=2)

par(mfrow=c(1,1),mar=c(10,4,3,1),oma=c(2,2,2,2))
boxplot(accuracy_df,names=c('class 1','class 2','total'),las=2,ylab='accuracy',ylim=c(0,1),
        main='Accuracy',col=colors[1:3])
legend ("topright", c(paste('regularization:',reg),paste('accuracy total:',round(mean(accuracy_df$accuracy)*100,2),'%')), lty=NULL)

boxplot(coef_variance_sorted_median,names=colnames(coef_variance_sorted_median),las=2,ylab='variebale importance sorted by median',
        main='Variable Importance with Regularization',col=colors[5])
legend ("topright", c(paste('regularization:',reg)), lty=NULL)

boxplot(coef_variance_sorted_mean,names=colnames(coef_variance_sorted_mean),las=2,ylab='variebale importance sorted by mean',
        main='Variable Importance with Regularization',col=colors[5])
legend ("topright", c(paste('regularization:',reg)), lty=NULL)

barplot(reduced_vector,las=2,ylab='mean of coefficients', main='Variable Importance with Regularization',col=colors[9])
legend ("topright", c(paste('regularization:',reg),paste('accuracy:',round(mean(accuracy_df$accuracy)*100,2),'%')), lty=NULL)

dev.off()

save(accuracy_df,file='/project/histone_marks_miRNA/logistic_regression/evaluation/pre_miRNA_promoter.RData')

save(accuracy_df,coef_variance_sorted_mean,file = paste(output_dir,reg,'_regularized_logistic_regression_training_test_accuracy_coefficients_',pdf_sub_name,'.RData',sep=''))

########################
#Prediction with CV:   #
########################

#parameters
runs = 20
alpha = 0.5 #1:Lasso, 0.5:ElasticNet, 0:Ridge
if(alpha == 0){reg = 'ridge'}; if(alpha == 0.5){reg = 'elnet'}; if(alpha == 1){reg = 'lasso'}
thr_coef = 0.05

#setup and run regularized logistic regression
coef_variance = as.data.frame(matrix(ncol=ncol(data_set),nrow=runs))
colnames(coef_variance) = colnames(data_set)
AUC = c()

for(run in 1:runs){
  print(run)
  best_lambdas = c()
  
  #run the model 200 times to get a sequence of 200 best lambdas
  for(i in 1:200){
    lmodel_binom.cv <- cv.glmnet(x = data.matrix(data_set), y = target, family='binomial',standardize=F,
                                 alpha=alpha,grouped=FALSE, nfolds=5, type.measure="auc")
    best_lambdas <- append(best_lambdas,lmodel_binom.cv$lambda.min)
  }
  
  #build model with sequence of best lambdas and get the best among them
  best_lambdas = sort(best_lambdas)
  lmodel_binom.cv <- cv.glmnet(x = data.matrix(data_set), y = target, family='binomial',lambda=best_lambdas,standardize=F,
                               alpha=alpha,grouped=FALSE, nfolds=5, type.measure="auc")
  
  best_lambda <- lmodel_binom.cv$lambda.min
  
  print(paste('optimal lambda:',round(log(best_lambda),2)))
  print(paste('auc:',round(max(lmodel_binom.cv$cvm),2)))
  
  #save coefficients and AUC
  coefs = coef(lmodel_binom.cv, s=best_lambda)
  coefs_vec = as.vector(coefs)
  coefs_vec = coefs_vec[2:length(coefs_vec)]
  coef_variance[run,] = coefs_vec
  
  AUC = c(AUC,round(max(lmodel_binom.cv$cvm),2))
  
}

#prepare coefficients for plotting
coef_variance_sorted = coef_variance

#shrinking outliers for plot
for(i in 1:ncol(coef_variance_sorted)){
  column = coef_variance_sorted[,i]
  b = boxplot(column)
  if(length(b$out)>0){
    for(j in 1:length(b$out)){
      if(b$out[j] > b$stats[5,1]){column[column==b$out[j]]=b$stats[5,1]}
      if(b$out[j] < b$stats[1,1]){column[column==b$out[j]]=b$stats[1,1]}
    }
  }
  coef_variance_sorted[,i] = column
}

#sort by decreasing mean
coef_variance_sorted = coef_variance_sorted[,c(names(sort(colMeans(coef_variance_sorted), decreasing=TRUE)))]

#if table with plotting names exist use them
if(exists('class_table')){
  old_names = names(coef_variance_sorted)
  new_names = c()
  for(name in old_names){
    new_name = class_table$abbr[class_table$abbr_GSE == name]
    if(length(new_name) == 0){
      new_name = name
    }
    new_names = append(new_names,new_name)
  }
  names(coef_variance_sorted) = new_names
  head(coef_variance_sorted,4)
}

#get mean coefficients and remove features with coefficients smaller than thr_coef
coefs_vec = colMeans(coef_variance_sorted)
names(coefs_vec) = colnames(coef_variance_sorted)
reduced_vector = coefs_vec[abs(coefs_vec) > thr_coef] 

#plot coefficients and AUC
estWidth <- min(4 * log(ncol(data_set)), 4 * log(nrow(data_set)))
estHeight <- min(4 * log(ncol(data_set)), 4 * log(nrow(data_set)))
CairoPDF(file = paste(output_dir,reg,'_regularized_logistic_regression_cv_',pdf_sub_name,'.pdf',sep=''), width = estWidth, height = estHeight)
par(mfrow=c(1,1),mar=c(10,4,3,1),oma=c(2,2,2,2))

boxplot(coef_variance_sorted,names=colnames(coef_variance_sorted),las=2,ylab='variebale importance sorted by mean',
        main='Variable Importance with Regularization',col=colors[5])
legend ("topright", c(paste('regularization:',reg),paste('AUC:',round(mean(AUC)*100,2),'%')), lty=NULL)

barplot(reduced_vector,las=2,ylab='mean of coefficients', main='Variable Importance with Regularization and 5 fold CV',col=colors[9])
legend ("topright", c(paste('regularization:',reg),paste('AUC:',round(mean(AUC)*100,2),'%')), lty=NULL)

dev.off() 
