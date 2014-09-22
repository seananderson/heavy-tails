# make a publication plot of p(nu < 20) by taxonomic order and class

# plot_nu_panel <- function

x <- filter(gomp_hat_ar1, taxonomic_class == "Mammalia") %>%
  arrange(taxonomic_class, order_p20) %>%


