# Weights-judging data reported in Davidson & Beaver (1977), p. 699

heaviness <- array(0, c(5,5,2), dimnames = list(
  c("90g","95g","100g","105g","110g"),
  c("90g","95g","100g","105g","110g"),
  c("order1", "order2") ))

# Order1: row stimulus first (in upper tri)
heaviness[,,1]  <- matrix(c( 0,  14,   6,   2,   1,
                            36,   0,  14,   7,   2,
                            44,  36,   0,  12,   5,
                            48,  43,  38,   0,  10,
                            49,  48,  45,  40,   0), 5,5, TRUE)

# Order2: column stimulus first (in upper tri)
heaviness[,,2]  <- matrix(c( 0,  18,  14,   3,   3,
                            32,   0,  16,   7,   4,
                            36,  34,   0,  16,  10,
                            47,  43,  34,   0,  22,
                            47,  46,  40,  28,   0), 5,5, TRUE)
