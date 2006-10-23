# Pork tasting data reported in Bradley & Terry (1952), p. 333

pork <- array(0, c(3,3,10), dimnames=list(c("C","Cp","CP"),
  c("C","Cp","CP"), c("jud1rep1","jud1rep2","jud1rep3",
  "jud1rep4","jud1rep5","jud2rep1","jud2rep2","jud2rep3",
  "jud2rep4","jud2rep5")))

# Judge 1, 5 repetitions
pork[,,1]  <- matrix(c(0,1,1,0,0,1,0,0,0), 3)
pork[,,2]  <- matrix(c(0,1,0,0,0,0,1,1,0), 3)
pork[,,3]  <- matrix(c(0,1,1,0,0,1,0,0,0), 3)
pork[,,4]  <- matrix(c(0,1,1,0,0,0,0,1,0), 3)
pork[,,5]  <- matrix(c(0,1,1,0,0,1,0,0,0), 3)

# Judge 2, 5 repetitions
pork[,,6]  <- matrix(c(0,1,0,0,0,0,1,1,0), 3)
pork[,,7]  <- matrix(c(0,1,0,0,0,0,1,1,0), 3)
pork[,,8]  <- matrix(c(0,0,0,1,0,1,1,0,0), 3)
pork[,,9]  <- matrix(c(0,0,0,1,0,1,1,0,0), 3)
pork[,,10] <- matrix(c(0,0,1,1,0,0,0,1,0), 3)
