s <- sprintf("s%X", 0:255)
while(length(s) > 1) {
	x <- seq(1,length(s),2)
	s <- sprintf("(%s:1,%s:1)", s[x],s[x+1])
}
s <- sprintf("%s;", s)
f <- file("balanced256.tree")
write(s,f)
close(f)


