aa <- "A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V"
aa <- strsplit(aa, "\\s+")[[1]]
m <- scan("wag.txt")
ss <- m[1:190]
p <- m[191:210]

s <- matrix(0,20,20)
s[upper.tri(s)] <- ss
s <- t(s)
s[upper.tri(s)] <- ss

o <- order(aa)
s <- s[o,][,o]

ss.new <- s[lower.tri(s)]

p.new <- p[o]

m <- matrix(sprintf("%.7f", ss.new),nrow=10)
m <- apply(m,2,paste, collapse=", ")
m <- paste(m, collapse=",\n")
cat(m); cat("\n\n")

m <- matrix(sprintf("%.7f", p.new),nrow=10)
m <- apply(m,2,paste, collapse=", ")
m <- paste(m, collapse=",\n")
cat(m); cat("\n\n")

# p2 <- p.new/sum(p.new)
# s2 <- t(p2*s)*p2
# s2 <- s2/sum(s2)
# s2 <- s2/p2
# r2 <- rowSums(s2)
# r3 <- max(r2)
# diag(s2) <- r3-r2
# s2 <- s2/rowSums(s2)
