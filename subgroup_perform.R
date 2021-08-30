print('This sub_group function was conducted to illustrate subgroup analysis')
print('Only the coefficient of main vairable will be shown')
print('')
print('***********    conducted by Yaxin Luo  2783025361@qq.com     *******************')
print('')
print('The subgroup should be redefine in advance, named by _group as postfix')
print('')
print('For example :')
print('if the analysis would be conducted in female/male subgroup, a new variable should be named as sex_group')
print('Note: all the variable in parameter group shoulb be presented in the parameter covar')


subgroup_perform <- function(data, main, covar, outcome, type, group){
  temp <- data %>% 
    select(all_of(main),all_of(covar),all_of(outcome),matches('_group')) %>% 
    na.omit()
  for (i in colnames(temp)){
    if (is.character(temp[,i])) temp[,i] <- as.factor(temp[,i]) 
  }
  subtable <- NULL
  for (i in group){
    k <- paste0(i,'_group')
    subinfo <- NULL
    ana <- temp %>% 
      select(all_of(main),all_of(covar),k,all_of(outcome)) %>% 
      select(-all_of(i))
    colnames(ana)[which(colnames(ana)==k)] <- 'subgroup'
    colnames(ana)[which(colnames(ana)==main)] <- 'index'
    if (type %in% c('or','rr')) {
      colnames(ana)[which(colnames(ana)==outcome)] <- 'outcome'
      if (type == 'rr') ts <- glm(outcome~.+index*subgroup,family = poisson(),data = ana) %>% summary()
      if (type == 'or') ts <- glm(outcome~.+index*subgroup,family = binomial(),data = ana) %>% summary()
    }
    if (type == 'cox') {
      ana$suroutcome <-  Surv(ana[,outcome[1]],ana[,outcome[2]])
      ana <- ana %>% select(-all_of(outcome))
      colnames(ana)[which(colnames(ana)=='suroutcome')] <- 'outcome'
      if (type == 'cox') ts <- coxph(outcome~.+index*subgroup,data = ana) %>% summary()
    }
    ts <- ts$coefficients
    pint <- ifelse(ts[nrow(ts),'Pr(>|z|)']<0.001,'<0.001',ts[nrow(ts),'Pr(>|z|)'] %>% round(3))
    l <- levels(ana$subgroup)
    for (j in levels(ana$subgroup)){
      ana1 <- ana %>% filter(subgroup==j) %>% select(-subgroup)
      if (type =='rr') ts <- glm(outcome~.,family = poisson(),data = ana1) %>% summary() 
      if (type =='or') ts <- glm(outcome~.,family = binomial(),data = ana1) %>% summary()
      if (type =='cox') {
        ts <- coxph(outcome~.,data = ana1) %>% summary()
        nevent <- ts$nevent
      }
      ts <- ts$coefficients
      subinfo <- rbind(subinfo,
                       cbind(variable = j,
                             n = nrow(ana1),
                             event = ifelse(type %in% c('or','rr'),sum(ana1$outcome),nevent),
                             rr = paste0(exp(ts[2,1]) %>% round(2),'(',
                                         exp(ts[2,1]-1.96*ts[2,2]) %>% round(2),'-',
                                         exp(ts[2,1]+1.96*ts[2,2]) %>% round(2),')'),
                             p = ifelse(ts[2,'Pr(>|z|)']<0.001,
                                        '<0.001',
                                        ts[2,'Pr(>|z|)'] %>% round(3))))
    }
    subinfo <- rbind(c(i,'','','','',''),cbind(subinfo,c(pint,'')))
    subtable <- rbind(subtable,
                      subinfo)
  }
  subtable
}

