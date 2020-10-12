plot_inference_report <-  function(){



}


plot_loss <- function(x){




}




get_baf <- function(vcf) {
  tb1 <- vcf@gt %>% as_tibble()
  colnames(tb1) <-  c("FORMAT", "tmp")
  tb1 <- tb1 %>%  tidyr::separate(tmp, into = c("GT", "PL", "DP", "AD"), sep = ":")
  tb2 <- tb1 %>% select(AD) %>% separate(AD, into = "ref", "alt") %>%  mutate(cov = ref + alt, chr = )
}
