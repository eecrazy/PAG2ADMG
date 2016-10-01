source('full_admg_learning.R')

#Generate pag1 from AAAI17 paper
pag1 <- matrix(rep(0,16), nrow = 4)
pag1[3,1] <- pag1[4,2] <- 1
pag1[1,3] <- pag1[2,4] <- pag1[4,3] <- pag1[3,4] <- 2
pag1 <- make_pag_from_amat(pag1)

#Convert pag1 to the full set of ADMGs using pag2admg
admg_list1 <- pag2admg(pag1)

#Generate pag2 from AAAI17 paper
pag2 <- matrix(rep(0,16), nrow = 4)
pag2[3,1] <- pag2[3,2] <- 1
pag2[4,3] <- 3
pag2[3,4] <- pag2[1,3] <- pag2[2,3] <- 2
pag2 <- make_pag_from_amat(pag2)

#Convert pag1 to the full set of ADMGs using pag2admg
admg_list2 <- pag2admg(pag2)

#If you would like to plot the pags or admg's the following lines have been commented out
# plot(pag1)
# plot_all_admgs(admg_list1)
# plot(pag2)
# plot_all_admgs(admg_list2)