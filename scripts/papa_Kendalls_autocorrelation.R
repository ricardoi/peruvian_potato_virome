data(Safariland)
s1 <- strength(Safariland, type="Barrat")
s2 <- strength(Safariland, type="Bascompte")
plot(s1, s2, log="x")
cor.test(s1, s2, type="ken")

# Pearson's product-moment correlation
# 
# data:  s1 and s2
# t = 2.429, df = 25, p-value = 0.02266
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.06829817 0.70064004
# sample estimates:
#       cor 
# 0.4369709 

# for lower level:
strength(t(Safariland))

closeness_w(t(Safariland))

# dat.mat

p1 <- strength(dat.mat, type="Barrat")
p2 <- strength(dat.mat, type="Bascompte")
plot(p1, p2, log="x")
cor.test(p1, p2, type="ken")
# Pearson's product-moment correlation
# 
# data:  s1 and s2
# t = -0.74137, df = 8, p-value = 0.4797
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.7615942  0.4475188
# sample estimates:
# cor 
# -0.2535503 
# for lower level:
strength(dat.mat)
distance_w(dat.mat) #require(tnet)
closeness_w(dat.mat)

#  Kendall's Tau correlation correlation test for betweeness centrality
vir.byalt <- datb.sp.table$`lower level`
vir.byalt[c(11,12)]
vir.b <- vir.byalt$betweenness
vir.wb<- vir.byalt$weighted.betweenness
plot(vir.b, vir.wb)
cor.test(vir.b, vir.wb,  method="kendall", conf.level = 0.95)

#  Kendall's Tau correlation correlation test for closeness centrality
pap.byhost <- datb.sp.table$`higher level`
pap.byhost[c(13,14)]
pap.c <- pap.byhost$closeness
pap.wc<- pap.byhost$weighted.closeness
plot(pap.c, pap.wc)
cor.test(pap.c, pap.wc)


# Betweenness by specie :: to look at 
ranks <- sapply(c("nodf", "binmatnest", "wine", "sort"), function(x) 
  nestedrank(t(dat.mat), method=x)[[2]])
cor(ranks) # high correlation between sort and other indicate that only abundance matters
