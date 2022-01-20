library(MuMIn)
library(purrr)


# Have to manually load RAM database
cod<-taxonomy %>% 
   filter(family=="Gadidae")
sole<-taxonomy %>% 
  filter(family=="Soleidae")

all_sci_names<-c(cod$scientificname,sole$scientificname)

stock_ids<-stock %>% 
  filter(scientificname %in% all_sci_names) %>% 
  select(stockid)

stock_id_clean<-timeseries_values_views %>% 
  filter(stockid %in% stock_ids$stockid) %>%
  select(stockid,year,TBbest,TCbest) %>% 
  drop_na() %>% 
  group_by(stockid) %>% 
  summarise(diff=max(year)-min(year)) %>% 
  filter(diff>20)


remove_vec=c(1,6,9,12,19,21,22,28,42,51,52,55)

named_remove<-unique(surplus$stockid)[-remove_vec]  

Fish_data<-timeseries_values_views %>% 
  filter(stockid %in% stock_id_clean$stockid) %>% 
  filter(stockid %in% named_remove)

surplus<-Fish_data %>% 
  group_by(stockid) %>% 
  mutate(f_biomass=lead(TBbest)) %>% 
  mutate(surplus=f_biomass-TBbest+TCbest) %>% 
  select(stockid,year,surplus,TBbest,TCbest) %>% 
  drop_na()


## Now let fit a fox model to surplus data

fox<-function(m,carry,biomass){
 out= -2.718*m*(biomass/carry)*log(biomass/carry)
return(out)
}

one_stock<-surplus %>% 
  filter(stockid==stock_id_clean$stockid[8])

guess_vec=c(max(one_stock$biomass)*0.37,max(one_stock$biomass))

one_stock_nls=nls(surplus~fox(m,carry,biomass),
                  data=one_stock,
                  start=list(m=guess_vec[1],carry=guess_vec[2]),trace=TRUE )



## Setup a purrr demonstration with all the nls
all_nls_fcn<-function(surplus_df){
  nls(surplus~fox(m,carry,biomass),
  data=surplus_df,
  start=list(m=max(surplus_df$biomass)*0.37,carry=max(surplus_df$biomass)))
}


fox_all<-surplus_final %>%
  group_by(stockid) %>% 
  nest() %>% 
  mutate(nls_model=map(data,~all_nls_fcn(.x))) %>% 
  mutate(AICc=map(nls_model,AICc))

fox_all$nls_model[[1]]

AICc(one_stock_nls)


### Create random surplus function

r_surplus<-function(surplus,biomass){
  avg_sur=mean(surplus)
  sigma=sqrt(sum((surplus-avg_sur)^2)/length(surplus))
  
  likeli=1/sigma*sqrt(2*pi)*exp(-(surplus-avg_sur)^2/(2*sigma^2))

AICc= -2*log(prod(likeli))+2+4/(length(surplus)-1-1)
return(AICc)
}

r_surplus(one_stock$surplus,one_stock$TBbest)


r_avg<-function(surplus){
  avg_sur=mean(surplus)
  
  rmse=sqrt(mean((avg_sur-surplus)^2))
  
  return(rmse)
}

r_mse<-surplus %>%
  group_by(stockid) %>% 
  nest() %>% 
  mutate(RMSE=map_dbl(data,~r_avg(.x$surplus)))

## this works to get the AICc
r_all<-surplus %>%
  group_by(stockid) %>% 
  nest() %>% 
  mutate(AICc=map(data,~r_surplus(.x$surplus,.x$biomass))) %>% 
  unnest(data)

### Try tidyverse way of getting suprlus




### Likelihood stuff because it worked
liklihood<-function(x,surplus,biomass){

  m=x[1]
  carry=x[2]
  sigma=x[3]
  
  like=log(1/(sigma*sqrt(2*pi)))-(surplus-fox(m=m,carry=carry,biomass = biomass))^2/(2*sigma^2)
  
  
  
  out=-sum(like)
  
  return(out)
}


test=optim(par = c(321160,868000,sd(one_stock$surplus)),fn=liklihood,surplus=one_stock$surplus,biomass=one_stock$TBbest)

optim_fcn<-function(data){
  hold=optim(par=c(max(data$TBbest)*0.37,max(data$TBbest),sd(data$surplus)),fn=liklihood,surplus=data$surplus,biomass=data$TBbest)
  
  #out=hold$value+6+2*3*4/(length(data$TBbest)-3-1)
  
  
  return(hold$par)
}



fox_mle<-surplus %>%
  group_by(stockid) %>% 
  nest() %>% 
  mutate(AICc=map(data,~optim_fcn(.x))) 



fit<-bbmle::mle2(liklihood,start=list(m=209637,carry=1087270,sigma=0.03),data=list(surplus=one_stock$surplus,biomass=one_stock$TBbest))


test_fox<-function(x,surplus,biomass){
  m=x[1]
  carry=x[2]
  sigma=x[3]
  
  
  likeli=1/sigma*sqrt(2*pi)*exp(-(surplus-fox(m,carry,biomass))^2/(2*sigma^2))
  
  return(likeli)
}
test_fox(c(255000,1089722,0.03),one_model$surplus,one_model$TBbest)

fox_all$nls_model[[1]]
