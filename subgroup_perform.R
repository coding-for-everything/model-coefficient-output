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


subgroup_perform <- function(data, main, ref = NA,covar, outcome, type, group){
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
      select(all_of(main),all_of(covar),all_of(k),all_of(outcome)) %>% 
      select(-all_of(i)) %>% na.omit()
    colnames(ana)[which(colnames(ana)==k)] <- 'subgroup'
    colnames(ana)[which(colnames(ana)==main)] <- 'index'
    if (!is.na(ref)){ ana$index <- relevel(ana$index,ref = ref)}
    if (type %in% c('or','rr')) {
      colnames(ana)[which(colnames(ana)==outcome)] <- 'outcome'
      if (type == 'rr') {
        t1 <- glm(outcome~.+index*subgroup,family = poisson(),data = ana)
        t2 <- glm(outcome~.,family = poisson(),data = ana) 
      }
      if (type == 'or') {
        t1 <- glm(outcome~.+index*subgroup,family = binomial(),data = ana) 
        t2 <- glm(outcome~.,family = binomial(),data = ana) 
      }
      ts <- anova(t1,t2,test = 'Chisq')
      #print(ts$`Pr(>|Chi|)` )
      pint <- ts$`Pr(>Chi)` %>% na.omit()
    }
    if (type == 'cox') {
      ana$suroutcome <-  Surv(ana[,outcome[1]],ana[,outcome[2]])
      ana <- ana %>% select(-all_of(outcome))
      colnames(ana)[which(colnames(ana)=='suroutcome')] <- 'outcome'
      t1 <- coxph(outcome~.+index*subgroup,data = ana)
      t2 <- coxph(outcome~.,data = ana)
      ts <- anova(t1,t2)
      pint <- ts$`P(>|Chi|)` %>% na.omit()
    }

    pint <- ifelse(pint<0.001,'<0.001',pint %>% round(3))
    l <- levels(ana$subgroup)
    for (j in levels(ana$subgroup)){
      ana1 <- ana %>% filter(subgroup==j) %>% select(-subgroup)
      if (type =='rr') ts <- glm(outcome~.,family = poisson(),data = ana1) %>% summary() 
      if (type =='or') ts <- glm(outcome~.,family = binomial(),data = ana1) %>% summary()
      if (type =='cox') {
        ts <- coxph(outcome~.,data = ana1) %>% summary()
        nevent <- ts$nevent
      }
      ts <- ts$coefficients %>% as.data.frame()
      if (length(levels(ana1$index))==2 | is.numeric(ana1$index)){
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
      
      if (length(levels(ana1$index))>2){
        xx <- length(levels(ana1$index))
        number <- data.frame(index = levels(ana1$index))
        if (type %in% c('or','rr')){
          number <- left_join(number,
                              ana1 %>% group_by(index) %>% summarise(n = n(),event = sum(outcome)),
                              'index')
        }
        if (type == 'cox'){
          temp1 <- data.frame(index = ana1$index,
                              outcome = ana1$outcome %>% as.character()) %>% 
            mutate(outcome = ifelse(str_detect(outcome,'\\+'),0,1))
          number <- left_join(number,
                              temp1 %>% group_by(index) %>% summarise(n = n(),event = sum(outcome)),
                              'index')
        }
        temp1 <- which(str_detect(rownames(ts),'index'))
        if (type %in% c('rr','or')){
          ts <- ts[,c(1,2,ncol(ts))]
        }
        else{ ts <- ts[1,3,col(ts)]}
        subinfo <- rbind(subinfo,
                         cbind(variable = c(j,rep('',xx-1)),
                               number,
                               rr = c('Reference',paste0(exp(ts[temp1,1]) %>% round(2),'(',
                                           exp(ts[temp1,1]-1.96*ts[temp1,2]) %>% round(2),'-',
                                           exp(ts[temp1,1]+1.96*ts[temp1,2]) %>% round(2),')')),
                               p = c('',ifelse(ts[temp1,'Pr(>|z|)']<0.001,
                                          '<0.001',
                                          ts[temp1,'Pr(>|z|)'] %>% round(3)))))
      }

    }
    for (kk in colnames(subinfo)){
      if (is.factor(subinfo[,kk])) subinfo[,kk] <- as.character(subinfo[,kk] )
    }
    if (length(levels(ana1$index))==2 | is.numeric(ana1$index)) subinfo <- rbind(c(i,'','','','',''),cbind(subinfo,c(pint,'')))
    if (length(levels(ana1$index))>2){
      subinfo <- rbind(c(i,'','','','','',''),cbind(subinfo,pint = c(pint,rep('',2*xx-1))))
    }
    subtable <- rbind(subtable,
                      subinfo)
  }
  subtable
}

