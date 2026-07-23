num_pares <- ncol(H) / 2
resultado500 <- matrix(0, nrow = nrow(H), ncol = num_pares)

for (i in 1:num_pares) {
  col_impar <- 2*i - 1
  col_par <- 2*i
  resultado500[, i] <- ifelse(H500[, col_impar] < 0 & H500[, col_par] > 0, 1, 0)
}


nombres <- c("0", "-0.2", "-0.4", "-0.6", "-0.8", "-1")

colnames(resultado500) <- nombres

colSums(resultado500)


colSums(resultado10)
colSums(resultado25)
colSums(resultado50)
colSums(resultado75)
colSums(resultado100)
colSums(resultado500)

