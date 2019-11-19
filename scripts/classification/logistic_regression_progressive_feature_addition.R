##############################
#Libraries                   #
##############################

library(glmnet)
library(Cairo)
library(plyr)
library(RColorBrewer)
library(matrixStats)
colors = brewer.pal(9,'Set1')

##############################
#READ DATA                   #
##############################

output_dir = '/project/histone_marks_miRNA/logistic_regression/plots/HELA/'
pdf_sub_name = 'pre_miRNA_promoter_intragenic'

load('/project/histone_marks_miRNA/data/RData_classification/HELA/encode/pre_miRNA_promoter_class.RData')

head(data_set,4)

selected_features = c('m_beta_value','p_H3K36me3','m_H3K79me2','p_H3K79me2','p_H3K27ac','p_beta_value','p_H3K27me3',
                      'm_H3K9ac','m_H3K9me3')
data_set_selected = data_set[colnames(data_set) %in% selected_features]

##########################################################
#FUNCTIONS                                               #
##########################################################

#run logistic regression on training and test set
run_logistic_regression <- function(target,data_set_cv){
  accuracy = c()
  distance = c()
  true_pos_rate = c()
  true_neg_rate = c()
  
  #run logistic regression for 100 times to get reliable prediction accuracy
  for(i in 1:100){
    
    target_index_0 = grep(0,target)
    target_index_1 = grep(1,target)
    training_index = sort(c(sample(x=target_index_0, size=sample_size_training_data[1], replace=F),
                            sample(x=target_index_1, size=sample_size_training_data[2], replace=F)))
    
    cv_data = data_set_cv[training_index,,drop=F]
    cv_label = target[training_index]
    
    validation_data = data_set_cv[-training_index,,drop=F]
    validation_label = target[-training_index]
    
    log_reg = glm(cv_label ~ .,cv_data,family="binomial")
    p = predict(log_reg,newdata = validation_data, type = "link") #predict target as 1: p > 0 --> above decision boundary Wlr^T * X + W0 = 0
    TP = sum(validation_label[p > 0] == 1) 
    FN = sum(validation_label[p < 0] == 1)
    TN = sum(validation_label[p < 0] == 0) 
    FP = sum(validation_label[p > 0] == 0)
    true_pos_rate = c(true_pos_rate,(TP / (TP + FN)))
    true_neg_rate = c(true_neg_rate,(TN / (TN + FP)))
    accuracy = c(accuracy,mean(c(true_pos_rate[i],true_neg_rate[i])))
    distance = c(distance,round(abs(true_pos_rate[i]-true_neg_rate[i]),2))
  }
  results_frame = data.frame(accuracy = accuracy, true_pos_rate = true_pos_rate, true_neg_rate = true_neg_rate)
  reduced_results = results_frame[distance < thr_class,]
  return(reduced_results)
}


##########################################################
#LOGISTIC REGRESSION - progressive feature addition      #
##########################################################

#specify size of training sample for both classes
sample_size_training_data = c(700,700) #class: 0/1
thr_class = 0.7

#initialize the evaluation measures
accuracy_table = rep(0,ncol(data_set_selected))
names(accuracy_table) = colnames(data_set_selected)
true_neg_rate_table = rep(0,ncol(data_set_selected))
names(true_neg_rate_table) = colnames(data_set_selected)
true_pos_rate_table = rep(0,ncol(data_set_selected))
names(true_pos_rate_table) = colnames(data_set_selected)

#get best one-feature model
for(c in 1:ncol(data_set_selected)){
  data_set_cv = data_set_selected[,c,drop=F]
  results = run_logistic_regression(target,data_set_cv)
  
  if(nrow(results) == 0){
    accuracy_table[c] = 0
    true_pos_rate_table[c] = 0
    true_neg_rate_table[c] = 0
  }else{
    accuracy_table[c] = mean(results$accuracy)
    true_pos_rate_table[c] = mean(results$true_pos_rate)
    true_neg_rate_table[c] = mean(results$true_neg_rate)
  }
}

print(accuracy_table)

feature_max_acc = names(which.max(accuracy_table))
print(feature_max_acc)
print(max(accuracy_table))

fixed_features = data.frame(data_set_selected[,colnames(data_set_selected) == feature_max_acc])
colnames(fixed_features) = feature_max_acc
remaining_features = data_set_selected[,-which.max(accuracy_table)]

increase_in_acc = rep(0,(ncol(data_set_selected)-1))
increase_in_acc[1] = max(accuracy_table)

change_in_true_pos_rate = rep(0,(ncol(data_set_selected)-1))
change_in_true_pos_rate[1] = true_pos_rate_table[which.max(accuracy_table)]

change_in_true_neg_rate = rep(0,(ncol(data_set_selected)-1))
change_in_true_neg_rate[1] = true_neg_rate_table[which.max(accuracy_table)]

#progressively add remaining features
for(features in 2:(ncol(data_set_selected)-1)){
  accuracy_table = rep(0,ncol(remaining_features))
  names(accuracy_table) = colnames(remaining_features)
  true_neg_rate_table = rep(0,ncol(remaining_features))
  names(true_neg_rate_table) = colnames(remaining_features)
  true_pos_rate_table = rep(0,ncol(remaining_features))
  names(true_pos_rate_table) = colnames(remaining_features)
  
  for(c in 1:ncol(remaining_features)){
    data_set_cv = cbind(fixed_features,remaining_features[,c])
    colnames(data_set_cv) = c(colnames(fixed_features),colnames(remaining_features)[c])
    results = run_logistic_regression(target,data_set_cv)
    if(nrow(results) == 0){
      accuracy_table[c] = 0
      true_pos_rate_table[c] = 0
      true_neg_rate_table[c] = 0
    }else{
      accuracy_table[c] = mean(results$accuracy)
      true_pos_rate_table[c] = mean(results$true_pos_rate)
      true_neg_rate_table[c] = mean(results$true_neg_rate)
    }
    
  }
  print(accuracy_table)
  feature_max_acc = names(which.max(accuracy_table))
  print(feature_max_acc)
  print(max(accuracy_table))
  
  fixed_features = cbind(fixed_features,remaining_features[,colnames(remaining_features) == feature_max_acc])
  colnames(fixed_features)[features] = feature_max_acc
  remaining_features = remaining_features[,-which.max(accuracy_table)]
  
  increase_in_acc[features] = max(accuracy_table)
  change_in_true_pos_rate[features] = true_pos_rate_table[which.max(accuracy_table)]
  change_in_true_neg_rate[features] = true_neg_rate_table[which.max(accuracy_table)]
}

names(increase_in_acc) = colnames(fixed_features)
names(change_in_true_neg_rate) = colnames(fixed_features)
names(change_in_true_pos_rate) = colnames(fixed_features)

print(increase_in_acc)

#plot results

estWidth <- min(4 * log(ncol(data_set_selected)), 4 * log(nrow(data_set_selected)))
estHeight <- min(4 * log(ncol(data_set_selected)), 4 * log(nrow(data_set_selected)))
CairoPDF(file = 'test.pdf',width = estWidth, height = estHeight)

par(mfrow=c(1,1),mar=c(8,4,2,2),oma=c(2,2,2,2))
plot(increase_in_acc,type = 'l',xaxt='n',ylab = 'increase in accuracy', xlab='feature',
     main='progressive feature addition',ylim=c(0.2,0.9),col=colors[1])
lines(change_in_true_neg_rate,col=colors[2])
lines(change_in_true_pos_rate,col=colors[3])
legend ("topright", c("total", "class1", "class2"), col=c(colors[1:3]), lty=1, lwd=1)
axis(1,at=1:length(increase_in_acc),label=names(increase_in_acc),las=2)

log_reg = glm(target ~ .,fixed_features[,1:6],family="binomial")
barplot(sort(log_reg$coefficients[-1],decreasing=T),las=2,main='coefficient weights')


dev.off()

##########################################################
#LOGISTIC REGRESSION - performance vs model complexity   #
##########################################################

sample_size_training_data = c(1000,800) #class: 0/1
thr_class = 0.7

columns = 1:ncol(data_set_selected)
n_model_features = 6

increase_in_acc = rep(0,n_model_features)
change_in_true_pos_rate = rep(0,n_model_features)
change_in_true_neg_rate = rep(0,n_model_features)
best_combination = rep('none',n_model_features)

for(f in 1:n_model_features){
  
  #get all combinations of features
  df = data.frame(columns)
  if(f > 1){
    for(i in 2:f){
      df = cbind(df,columns)
    }
    combination = expand.grid(df)
    
    duplicated_rows_to_remove = c()
    for(r in 1:nrow(combination)){
      row = combination[r,]
      if(sum(duplicated(unlist(row))) > 0){
        duplicated_rows_to_remove = c(duplicated_rows_to_remove,r)
      }
    }
    combination = combination[-duplicated_rows_to_remove,]
    
    for(r in 1:nrow(combination)){
      row = combination[r,]
      combination[r,] = sort(row)
    }
    combination = unique(combination)
    rownames(combination) = 1:nrow(combination)
  }else{
    combination = df
  }
  
  coln = combination[,1]
  if(f > 1){
    for(c in 2:ncol(combination)){
      coln = paste(coln,combination[,c])
    }
  }

  accuracy_table = rep(0,nrow(combination))
  names(accuracy_table) = coln
  true_neg_rate_table = rep(0,nrow(combination))
  names(true_neg_rate_table) = coln
  true_pos_rate_table = rep(0,nrow(combination))
  names(true_pos_rate_table) = coln
  
  for(comb in 1:nrow(combination)){
    feature_combination = unlist(combination[comb,])
    data_set_cv = data_set_selected[,feature_combination,drop=F]
    results = run_logistic_regression(target,data_set_cv)
    if(nrow(results) == 0){
      accuracy_table[comb] = 0
      true_pos_rate_table[comb] = 0
      true_neg_rate_table[comb] = 0
    }else{
      accuracy_table[comb] = mean(results$accuracy)
      true_pos_rate_table[comb] = mean(results$true_pos_rate)
      true_neg_rate_table[comb] = mean(results$true_neg_rate)
    }
  }
  print(accuracy_table)
  feature_max_acc = names(which.max(accuracy_table))
  print(feature_max_acc)
  print(max(accuracy_table))
  
  increase_in_acc[f] = max(accuracy_table)
  change_in_true_pos_rate[f] = true_pos_rate_table[which.max(accuracy_table)]
  change_in_true_neg_rate[f] = true_neg_rate_table[which.max(accuracy_table)]
  best_combination[f] = feature_max_acc
}


estWidth <- min(4 * log(ncol(data_set_selected)), 4 * log(nrow(data_set_selected)))
estHeight <- min(4 * log(ncol(data_set_selected)), 4 * log(nrow(data_set_selected)))
CairoPDF(file = paste(output_dir,'accuracy_vs_model_complexity_',pdf_sub_name,'.pdf',sep=''),width = estWidth, height = estHeight)

par(mfrow=c(1,1),mar=c(8,4,2,2),oma=c(2,2,2,2))
plot(increase_in_acc,type = 'l',xaxt='n',ylab = 'increase in accuracy', xlab='feature',
     main='progressive feature addition',ylim=c(0.2,0.9),col=colors[1])
lines(change_in_true_neg_rate,col=colors[2])
lines(change_in_true_pos_rate,col=colors[3])
legend.text = c("total", "class1", "class2",colnames(data_set_selected))
legend ("bottomright", legend.text, col=c(colors[1:3],rep('black',ncol(data_set_selected))), pch=c('*','*','*',1:ncol(data_set_selected)),ncol=3)
axis(1,at=1:length(increase_in_acc),label=best_combination,las=2)

dev.off()

##########################################################
#REGULARIZED LOGISTIC REGRESSION - on selected features  #
##########################################################

#run a regularized logistic regression on selected features
best_features = c(as.numeric(unlist(strsplit(best_combination[5],' '))),2,4)#-> to include both features for methylation
data_set = data_set_selected[,best_features]

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
  
  lmodel_binom = glmnet(x = data.matrix(training_data_set), y = training_target, family='binomial',standardize=F,lambda=best_lambdas,alpha=alpha) 
  predictions = predict(lmodel_binom, data.matrix(test_data_set), type="class")
  
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



