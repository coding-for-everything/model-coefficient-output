full_glm <- function(data,main,ref = NULL,other = NULL,merger = FALSE,covar = NULL,outcome,type){
  if (is.null(covar)) covar <- main
  data <- data %>% 
    dplyr::select(all_of(main),all_of(covar),all_of(outcome))
  colnames(data)[which(colnames(data)==main)] <- 'index'
  colnames(data)[which(colnames(data)==outcome)] <- 'outcome'
  data <- na.omit(data)
  
  temp <- data.frame(n = nrow(data),
                     event = sum(data$outcome)) %>% as.data.frame()
  
  if (!is.numeric(data$index)){
    if (is.null(ref)|is.null(other)) return('Please clarify parameters: ref/other')
    
    data <- data %>% 
      filter(index %in% c(all_of(ref),all_of(other))) %>% 
      mutate(index = as.factor(index) %>% relevel(ref = all_of(ref)))
    # if (merge == TRUE) data$index <- ifelse(index == ref,ref,other[1]) %>% as.factor() %>% relevel(ref = all_of(ref))
    
    temp <- data %>% 
      group_by(index) %>% 
      summarise(n = n(),
                event = ifelse(is.numeric(outcome),sum(outcome),sum(as.numeric(str_extract(outcome,'[0-9]')))))# %>% 
    # cbind(table(data$index %>% as.character()),.) %>% as.data.frame()
  }
  
  # colnames(temp)[1] <- 'index'
  for (i in colnames(data)){
    if (is.character(data[,i])) data[,i] <- as.factor(data[,i] %>% unlist())
  }
  
  if (type == 'rr') t <- glm(outcome~.,data = data,family = 'poisson') %>% summary()
  if (type == 'or') t <- glm(outcome~.,data = data,family = 'binomial') %>% summary()
  
  t <- t$coefficients %>% as.data.frame() %>% .[-1,]
  table <- cbind(risk = paste0(exp(t[,1]) %>% round(2),' (',
                               exp(t[,1] - 1.96*t[,2]) %>% round(2),'-',
                               exp(t[,1] + 1.96*t[,2]) %>% round(2),')'),
                 pvalue = ifelse((t[,'Pr(>|z|)'] %>% round(3))==0,'<0.001',t[,'Pr(>|z|)'] %>% round(3)))
  
  showcase <- NULL
  j <- 0
  for (i in colnames(data %>% dplyr::select(-outcome))){
    j <- j+1
    showcase <- rbind(showcase,c(i,'',''))
    if (is.numeric(data[,i])){
      showcase <- rbind(showcase,c(paste0('   ',i),table[j,]))
      next()
    }
    showcase <- rbind(showcase,
                      cbind(levels(data[,i]) %>% unlist() %>% paste0('   ',.),
                            c('Reference',table[,'risk'][j:(j+length(levels(data[,i]))-2)]),
                            c('',table[,'pvalue'][j:(j+length(levels(data[,i]))-2)])))
    j <- j+length(levels(data[,i]))-2
  }
  colnames(showcase) <- c('Variable',type,'pvalue')
  showcase[1,1] <- main
  cbind(showcase,
        n = c('',temp$n,rep('',nrow(showcase)-nrow(temp)-1)),
        event = c('',temp$event,rep('',nrow(showcase)-nrow(temp)-1))) 
}

full_cox <- function(stime,covar,main,data){
  data <- data %>% 
    dplyr::select(all_of(stime),all_of(main),all_of(covar)) %>% 
    na.omit()
  for (i in 1:ncol(data)){
    if (is.character(data[,i]%>% unlist()) ) data[,i] <- as.factor(data[,i] %>% unlist())
  }
  s <- Surv(data[,stime[1]] %>% unlist(),data[,stime[2]] %>% unlist)
  data <- data %>% dplyr::select(-c(all_of(stime)))
  sx <- coxph(s~.,data = data) %>% summary()
  ss <- sx$conf.int
  hrtable <- cbind(hr = paste0(ss[,1] %>% round(2),'(',ss[,3] %>% round(2),'-',ss[,4] %>% round(2),')'),
                   pvalue = ifelse(sx$coefficients[,'Pr(>|z|)'] %>% round(3) == 0,'<0.001',sx$coefficients[,'Pr(>|z|)'] %>% round(3)))
  showcase <- NULL
  j <- 0
  for (i in colnames(data)){
    j <- j+1
    showcase <- rbind(showcase,c(i,'',''))
    if (is.numeric(data[,i] %>% unlist())) {
      showcase <- rbind(showcase,c(i,hrtable[j,]))
      next()
    }
    showcase <- rbind(showcase,
                      cbind(levels(data[,i] %>% unlist()),
                            c('Reference',hrtable[,'hr'][j:(j+length(levels(data[,i] %>% unlist()))-2)]),
                            c('',hrtable[,'pvalue'][j:(j+length(levels(data[,i] %>% unlist()))-2)])))
    j <- j+length(levels(data[,i] %>% unlist()))-2
  } 
  cbind(showcase,
        concordance = c(sx$concordance[1] %>% round(3),rep('',nrow(showcase)-1))) 
}

uniq_glm <- function(data, main, covar = NULL, outcome, type){
  data <- data %>% select(all_of(main),all_of(covar),all_of(outcome))
  if (is.null(covar)) covar <- main
  for (i in c(main,covar)) {
    if (!is.numeric(data[,i] %>% unlist())) data[,i] <- as.factor(data[,i])
  }
  showcase <- data.frame()
  for (i in c(main,covar)){
    temp <- data %>% dplyr::select(all_of(i),all_of(outcome))
    colnames(temp) <- c('index','outcome')
    
    if (is.factor(temp$index)){
      num <- temp %>% na.omit() %>%  group_by(index) %>% 
        summarise(n = n(),
                  nevent = sum(outcome))
    } else {
      num <- data.frame(n = nrow(na.omit(temp)),
                        nevent = sum(na.omit(temp)[,'outcome']))
    } 
    
    if (type =='rr') t <- glm(outcome~index,data = temp,family = 'poisson') %>% summary()
    if (type == 'or') t <- glm(outcome~index,data = temp,family = 'binomial') %>% summary()
    
    t <- t$coefficients %>% as.data.frame() %>% .[-1,]
    table <- cbind(risk = paste0(exp(t[,1]) %>% round(2),' (',
                                 exp(t[,1] - 1.96*t[,2]) %>% round(2),'-',
                                 exp(t[,1] + 1.96*t[,2]) %>% round(2),')'),
                   pvalue = ifelse((t[,'Pr(>|z|)'] %>% round(3))==0,'<0.001',t[,'Pr(>|z|)'] %>% round(3)))
    showcase <- rbind(showcase,data.frame(variable = all_of(i),
                                          risk = '',
                                          pvalue = '',
                                          n = '',
                                          nevent = '',stringsAsFactors = FALSE),stringsAsFactors = FALSE)
    if (is.numeric(temp$index)) {
      table <- cbind(variable = paste0('  ',all_of(i)) %>% as.character(),
                     table ,
                     num,stringsAsFactors = FALSE) 
    }
    if (is.factor(temp$index)){
      table <- rbind(c('Reference',''),table)
      table <- cbind(table,num[,-1],stringsAsFactors = FALSE)
      table <- cbind(levels(temp$index),table)
      colnames(table) <- colnames(showcase)
    }
    showcase <- rbind(showcase,table,stringsAsFactors = FALSE)
  }
  showcase
} 

uniq_cox <- function(data, stime, outcome, main, covar = NULL){
  if (is.null(covar)) covar <- all_of(main)
  for (i in c(main,covar)) {
    if (!is.numeric(data[,i] %>% unlist())) data[,i] <- as.factor(data[,i])
  }
  
  showcase <- data.frame()
  for (i in c(main,covar)){
    temp <- data %>% dplyr::select(all_of(i),all_of(stime),all_of(outcome))
    colnames(temp) <- c('index','stime','outcome')
    
    if (is.factor(temp$index)){
      num <- temp %>% na.omit() %>% 
        group_by(index) %>% 
        summarise(n = n(),
                  nevent = sum(outcome))
    } else {
      num <- data.frame(n = nrow(na.omit(temp)),
                        nevent = sum(na.omit(temp)[,'outcome']))
    } 
    t <- coxph(Surv(stime,outcome)~index,data = temp) %>% summary()
    t <- t$coefficients %>% as.data.frame() 
    table <- cbind(hr = paste0(exp(t[,1]) %>% round(2),' (',
                               exp(t[,1] - 1.96*t[,2]) %>% round(2),'-',
                               exp(t[,1] + 1.96*t[,2]) %>% round(2),')'),
                   pvalue = ifelse((t[,'Pr(>|z|)'] %>% round(3))==0,'<0.001',t[,'Pr(>|z|)'] %>% round(3)))
    showcase <- rbind(showcase,data.frame(variable = all_of(i),
                                          hr = '',
                                          pvalue = '',
                                          n = '',
                                          nevent = '',stringsAsFactors = FALSE),stringsAsFactors = FALSE)
    if (is.numeric(temp$index)) {
      table <- cbind(variable = paste0('  ',all_of(i)) %>% as.character(),
                     table ,
                     num,stringsAsFactors = FALSE) 
    }
    if (is.factor(temp$index)){
      table <- rbind(c('Reference',''),table)
      table <- cbind(table,num[,-1],stringsAsFactors = FALSE)
      table <- cbind(levels(temp$index),table)
      colnames(table) <- colnames(showcase)
    }
    showcase <- rbind(showcase,table,stringsAsFactors = FALSE)
  }
  showcase
}
print('There are four function for glm regression and cox regression, ')
print('Two for univairable analysis, named as uniq_glm and uniq_cox') 
print('Two for multivariable analysis,named as full_glm  and  ful_cox')
print('           *****************           ')
print('Conducted by Yaxin Luo, 2783025361@qq.com')
print('           *****************           ')
print('')
print('')
print('For univariable analysis: ')
print('parameter includes data, main, covar for both')
print('logistic includes parameters outcome, while cox includes parameters stime and outcome')
print('')
print('Note: as the function name shows, these uniq_ functions were used for univariable analysis')
print('Showing in one sheet just for convinence, for multivariable analysis, please use full_ function')
