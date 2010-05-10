## Perceived health risk of drugs (Wickelmaier, useR! 2008)

drugrisk <- array(0, c(6,6,4), dimnames = list(
  c("alc", "tob", "can", "ecs", "her", "coc"),
  c("alc", "tob", "can", "ecs", "her", "coc"),
  c("female30", "female31", "male30", "male31") ))
names(dimnames(drugrisk)) <- c(">", "<", "group")

## Group1: female participants, age <= 30
drugrisk[,,1]  <- matrix(c( 0, 25, 32,  6,  3,  4,
                           23,  0, 29,  2,  2,  3,
                           16, 19,  0,  0,  2,  3,
                           42, 46, 48,  0,  8, 22,
                           45, 46, 46, 40,  0, 38,
                           44, 45, 45, 26, 10,  0), 6, 6, TRUE)

## Group2: female participants, age >= 31
drugrisk[,,2]  <- matrix(c( 0, 32, 27, 12,  0,  8,
                           16,  0, 14,  5,  2,  6,
                           21, 34,  0, 10,  3,  4,
                           36, 43, 38,  0,  8, 20,
                           48, 46, 45, 40,  0, 42,
                           41, 42, 44, 28,  6,  0), 6, 6, TRUE)

## Group3: male participants, age <= 30
drugrisk[,,3]  <- matrix(c( 0, 28, 35, 10,  4,  7,
                           20,  0, 18,  2,  0,  3,
                           13, 30,  0,  3,  1,  0,
                           38, 46, 45,  0,  1, 17,
                           44, 48, 47, 47,  0, 44,
                           41, 45, 48, 31,  4,  0), 6, 6, TRUE)

## Group4: male participants, age >= 31
drugrisk[,,4]  <- matrix(c( 0, 27, 30, 10,  4,  9,
                           21,  0, 21,  7,  6,  9,
                           18, 27,  0,  6,  2,  2,
                           38, 41, 42,  0,  4, 30,
                           44, 42, 46, 44,  0, 41,
                           39, 39, 46, 18,  7,  0), 6, 6, TRUE)
