## code to prepare `DATASET` dataset goes here

set.seed(100)
HomUHet_data=HomUHet.sim(Pred_type="Gaussian", J=500, K=4, beta=NULL,
                      rho=0.5,sigma=2, level="m",
                      nlower=50,nupper=300)
usethis::use_data(HomUHet_data, overwrite = TRUE)
