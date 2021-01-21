library(tidyverse)

discretize <- function(v,t, p_val= 10^-6,alpha_fp=0.001,alpha_fn=0.001){
  if(t==0){
    return(3)
  }
  b.test <- binom.test(v,t,alpha_fp,alternative="g")
  if(b.test$p.value > p_val){
    return(0)
  }
  b.test <- binom.test(v,t,1-alpha_fn, alternative="l")
  if(b.test$p.value > p_val){
    return(2)
  }else{
    return(1)
  }
}



ad <- snakemake@input[['ad_file']]
dp <- snakemake@input[['dp_file']]
pval <- snakemake@params[['pval']]
alpha_fp <- snakemake@params[['alpha_fp']]
alpha_fn <- snakemake@params[['alpha_fn']]


#test code 
# pth <- "/scratch/data/leah/doubletD/test"
# ad <- file.path(pth, "AD.csv")
# dp <-file.path(pth, "DP.csv")
# 
# pval <- 10^-6
# alpha_fp <- 0.001
# alpha_fn <- 0.001

ad.dat <- read.csv(ad) %>% pivot_longer(cols=c(contains("var"), contains("chr")), names_to="variant", values_to="ad")
dp.dat <- read.csv(dp)%>% pivot_longer(cols=c(contains("var"), contains("chr")), names_to="variant", values_to="dp")

dat <- inner_join(ad.dat, dp.dat, by=c("cell_id", "variant")) %>% 
  rowwise() %>%
  mutate(zygosity= discretize(ad, dp, pval, alpha_fp, alpha_fn))%>%
  ungroup() %>%
  select(cell_id, variant, zygosity) %>%
  pivot_wider(names_from=variant, values_from=zygosity)

write.table(dat, snakemake@output[['outfile']], sep="\t", row.names=F)








# discretize(7,0)
# b.test <- binom.test(75, 100,1-0.001,alternative="l")
