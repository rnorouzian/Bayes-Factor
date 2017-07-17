
D = read.csv("C:/Users/Reza/Google Drive/Manuscript 2/BFdPvalue.csv")


BF.d.pvalue = Vectorize(function(t, n1, n2 = NA, scale = sqrt(2)/2, log.BF = FALSE){
  
    options(warn = -1)  
      t = abs(t)
      N = ifelse(is.na(n2), n1, (n1*n2)/(n1+n2))
     df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
      d = t / sqrt(N)
  
     H1 = integrate(function(delta) dcauchy(delta, 0, scale) * dt(t, df, delta*sqrt(N)), -Inf, Inf)[[1]]
     H0 = dt(t, df)
   BF10 = ifelse(log.BF, log(H1/H0), H1/H0)
p.value = 2*(1-pt(t, df))
  
  cbind(BF10 = BF10, p.value = p.value, d = d, H1 = H1, H0 = H0)
  
}, vectorize.args = c("t", "n1", "n2", "scale", "log.BF"))

  
  
t.value = D$t.value ; n1 = D$n1 ; n2 = D$n2 ; t.type = D$t.type ; p.value = D$p.value ; Data.size = nrow(D)


## Pre-analysis: 
pre_analysis = function(){

list(Data.size = Data.size,
paired.samples = sum(t.type == 0), #Or length(t.type[t.type == 0]) table(t.type <= 1)[c("TRUE", "FALSE")]
independent.samples = sum(t.type == 1),
sig.pvalues = unname(table(p.value <= .05)["TRUE"]),
non.sig.pvalues = unname(table(p.value <= .05)["FALSE"]))
  
}



simulation = function(){
  
t.sim = mapply(FUN = rt, n = 1, df = 2:1e3, ncp = seq(0, 1, l = 999))
b = BF.d.pvalue(t = t.sim, n1 = 3:1001, n2 = NA)

pairs(b[2 , ]~ b[3 , ] + b[1 , ], data = b, gap = .2, t = "n", diag.panel = panel.hist,
      lower.panel = NULL, las = 1, panel = panel.smooth, 
      labels = c("p-value", "Cohen's d" , "BF10"), col = 3)
}



## Main Analysis:
b = BF.d.pvalue(t = t.value, n1 = n1, n2 = n2)

BF = b[1, ]  ;  p.value = b[2, ]   ;  d = b[3, ]


pval.bf = data.frame(p = p.value, BF = BF) 
d.bf = data.frame(BF = BF, d = d)  
pval.d = data.frame(p = p.value, d = d)



p.value_BF = function(){
  
## Subsetting for a 7x4 plot p.value against BF:
p.b7.4 = subset(pval.bf, (BF > .1 & BF <= 1/3) & (p > .05)) 
p.b6.4 = subset(pval.bf, (BF > 1/3 & BF <= 1) & (p > .05))
p.b5.4 = subset(pval.bf, (BF > 1 & BF <= 3) & (p > .05)) 
p.b4.4 = subset(pval.bf, (BF > 3 & BF <= 10) & (p > .05)) 
p.b3.4 = subset(pval.bf, (BF > 10 & BF <= 30) & (p > .05))  
p.b2.4 = subset(pval.bf, (BF > 30 & BF <= 100) & (p > .05)) 
p.b1.4 = subset(pval.bf, (BF > 100) & (p > .05)) 


p.b7.3 = subset(pval.bf, (BF > .1 & BF <= 1/3) & (p > .01 & p <= .05)) 
p.b6.3 = subset(pval.bf, (BF > 1/3 & BF <= 1) & (p > .01 & p <= .05))
p.b5.3 = subset(pval.bf, (BF > 1 & BF <= 3) & (p > .01 & p <= .05)) 
p.b4.3 = subset(pval.bf, (BF > 3 & BF <= 10) & (p > .01 & p <= .05))
p.b3.3 = subset(pval.bf, (BF > 10 & BF <= 30) & (p > .01 & p <= .05))
p.b2.3 = subset(pval.bf, (BF > 30 & BF <= 100) & (p > .01 & p <= .05))
p.b1.3 = subset(pval.bf, (BF > 100) & (p > .01 & p <= .05))


p.b7.2 = subset(pval.bf, (BF > .1 & BF <= 1/3) & (p > .001 & p <= .01))
p.b6.2 = subset(pval.bf, (BF > 1/3 & BF <= 1) & (p > .001 & p <= .01))
p.b5.2 = subset(pval.bf, (BF > 1 & BF <= 3) & (p > .001 & p <= .01))
p.b4.2 = subset(pval.bf, (BF > 3 & BF <= 10) & (p > .001 & p <= .01))
p.b3.2 = subset(pval.bf, (BF > 10 & BF <= 30) & (p > .001 & p <= .01))
p.b2.2 = subset(pval.bf, (BF > 30 & BF <= 100) & (p > .001 & p <= .01))
p.b1.2 = subset(pval.bf, (BF > 100) & (p > .001 & p <= .01))


p.b7.1 = subset(pval.bf, (BF > .1 & BF <= 1/3) & (p <= .001))
p.b6.1 = subset(pval.bf, (BF > 1/3 & BF <= 1) & (p <= .001))
p.b5.1 = subset(pval.bf, (BF > 1 & BF <= 3) & (p <= .001))
p.b4.1 = subset(pval.bf, (BF > 3 & BF <= 10) & (p <= .001))
p.b3.1 = subset(pval.bf, (BF > 10 & BF <= 30) & (p <= .001))
p.b2.1 = subset(pval.bf, (BF > 30 & BF <= 100) & (p <= .001))
p.b1.1 = subset(pval.bf, (BF > 100) & (p <= .001))


pb.cat.sizes = lapply(mget(ls(pattern = "p\\.b\\d")), nrow)
total.pb.cat.sizes = Data.size


bf.Decisive = nrow(subset(pval.bf, (BF > 100)))
bf.very.strong = nrow(subset(pval.bf, (BF > 30 & BF <= 100)))
bf.strong = nrow(subset(pval.bf, (BF > 10 & BF <= 30)))
bf.substan.H1 = nrow(subset(pval.bf, (BF > 3 & BF <= 10)))
bf.Anecd.H1 = nrow(subset(pval.bf, (BF > 1 & BF <= 3)))
bf.Anecd.H0 = nrow(subset(pval.bf, (BF > 1/3 & BF <= 1)))
bf.substan.H0 = nrow(subset(pval.bf, (BF > .1 & BF <= 1/3)))

bf.cat.check = sum(bf.Decisive, bf.very.strong, bf.strong, bf.substan.H1, bf.Anecd.H1, bf.Anecd.H0, bf.substan.H0)


p.Decisive = nrow(subset(pval.bf, p <= .001)) 
p.substan = nrow(subset(pval.bf, (p > .001 & p <= .01)))
p.positive = nrow(subset(pval.bf, (p > .01 & p <= .05)))
p.no = nrow(subset(pval.bf, p > .05))

p.cat.check = sum(p.Decisive, p.substan, p.positive, p.no)
  
  
result = matrix(unlist(pb.cat.sizes), ncol = 4, byrow = TRUE)
dimnames(result) = list(BF = c("Decisive", "Very Strong", "Strong", "Substantial", "Anecdotal", "Anecdotal(H0)", "Substantial(H0)")
                             ,p.value = c("Decisive", "Substantial", "Positive", "None"))

p.marginal = matrix(c(p.Decisive, p.substan, p.positive, p.no), nrow = 1)
bf.marginal = matrix(c(bf.Decisive, bf.very.strong, bf.strong, bf.substan.H1,
                       bf.Anecd.H1, bf.Anecd.H0, bf.substan.H0, bf.cat.check), ncol = 1)

X.sq = chisq.test(result, correct = FALSE)

chisq.table = cbind(rbind(result, p.marginal), bf.marginal)
dimnames(chisq.table) = list(BF = c("Decisive", "Very Strong", "Strong", "Substantial", "Anecdotal", "Anecdotal(H0)", "Substantial(H0)", "Marginal.p")
                             ,p.value = c("Decisive", "Substantial", "Positive", "None", "Marginal.BF"))


# Regular expression to create all 28 plot names:
plot.names = noquote(sprintf("p.b%d.%d", 1:7, rep(1:4, each = 7)))


# Split screen:
original_par = par(no.readonly = TRUE)
on.exit(par(original_par))

png("Plot_pval_BF.png", res = 500, width = 4.5, height = 5.3, units = "in")
par(mai = rep(.5, 4)) 
par(mfcol = c(7, 4), mar = rep(.07, 4), oma = c(4.5, 5, 4.2, 5.1))


#for graphic device only use: par(mfcol = c(7, 4), mar = rep(.08, 4), oma = rep(7, 4))

## Multiple plotting:
invisible(lapply(mget(plot.names), function(pb) 
  if(nrow(pb) == 0){
    plot.new(); box()
  }else{
    plot(pb, pch = 21, bg = 3, cex = 1.3, xaxt = "n", yaxt = "n")
  }))



## Advanced Labeling:

mid.x = seq(grconvertX(0 + (1 /  8), "nic"), grconvertX(1 - (1 /  8), "nic"), l = 4)
mid.y = seq(grconvertY(0 + (1 / 14), "nic"), grconvertY(1 - (1 / 14), "nic"), l = 7)

gap.x = sort(c(mid.x[-length(mid.x)] + diff(mid.x)[1L] / 2, grconvertX(0:1, "nic")))
gap.y = sort(c(mid.y[-length(mid.y)] + diff(mid.y)[1L] / 2, grconvertY(0:1, "nic")))

xlim = c( grconvertX(0, "nic"), grconvertX(1, "nic") )
ylim = c( grconvertY(0, "nic"), grconvertY(1, "nic") )


v = c(Substantial = bf.substan.H0, Anecdotal = bf.Anecd.H0, Anecdotal = bf.Anecd.H1, Substantial = bf.substan.H1, 
      Strong=bf.strong, "Very Strong" = bf.very.strong, Decisive = bf.Decisive)/Data.size*1e2
l = paste0(names(v), "\n", v, "%")       

mtext(l, side = 4, at = mid.y, las = 1, line = .2, cex = .7, font = 2)


v = c(Decisive= p.Decisive, Substantial= p.substan, Positive= p.positive, None= p.no)/Data.size*1e2
l = paste0(names(v), "\n", v, "%")

text(mid.x, ylim[2] + .05, l, xpd = NA, cex = 1.1, font = 2)


l = c('0', '.001', '.01', '.05', '1') 
mtext(l, side = 1, at = gap.x, line = .55, cex = .9, font = 2)


l = c('1/10', '1/3', '1', '3', '10', '30', '100', '> 100000')
text(xlim[1], gap.y, l, xpd = NA, adj = 1.05, cex = 1.3, font = 2)

mtext(expression(bolditalic("p")*bold("-value")), side = 1, at = mean(gap.x), line = 3.5, cex = 1.2)
text(xlim[1], mean(gap.y), "Bayes Factor", xpd = NA, adj = c(.5, -2.2), cex = 1.8, font = 2, srt = 90)
text(mean(gap.x), ylim[2] + .16, bquote(bold("Evidence against"~bolditalic(H)[0])), cex = 1.8, xpd = NA)
 
dev.off()

list(chisq.table = chisq.table, X.sq = X.sq)

}




d_BF = function(){
  
  ## Subsetting for a 7x4 plot Cohen's d against BF:
  d.b7.4 = subset(d.bf, (BF > .1 & BF <= 1/3) & (d > .8)) 
  d.b6.4 = subset(d.bf, (BF > 1/3 & BF <= 1) & (d > .8))
  d.b5.4 = subset(d.bf, (BF > 1 & BF <= 3) & (d > .8)) 
  d.b4.4 = subset(d.bf, (BF > 3 & BF <= 10) & (d > .8)) 
  d.b3.4 = subset(d.bf, (BF > 10 & BF <= 30) & (d > .8))  
  d.b2.4 = subset(d.bf, (BF > 30 & BF <= 100) & (d > .8)) 
  d.b1.4 = subset(d.bf, (BF > 100) & (d > .8)) 
  
  
  d.b7.3 = subset(d.bf, (BF > .1 & BF <= 1/3) & (d > .5 & d <= .8)) 
  d.b6.3 = subset(d.bf, (BF > 1/3 & BF <= 1) & (d > .5 & d <= .8))
  d.b5.3 = subset(d.bf, (BF > 1 & BF <= 3) & (d > .5 & d <= .8)) 
  d.b4.3 = subset(d.bf, (BF > 3 & BF <= 10) & (d > .5 & d <= .8))
  d.b3.3 = subset(d.bf, (BF > 10 & BF <= 30) & (d > .5 & d <= .8))
  d.b2.3 = subset(d.bf, (BF > 30 & BF <= 100) & (d > .5 & d <= .8))
  d.b1.3 = subset(d.bf, (BF > 100) & (d > .5 & d <= .8))
  
  
  d.b7.2 = subset(d.bf, (BF > .1 & BF <= 1/3) & (d > .2 & d <= .5))
  d.b6.2 = subset(d.bf, (BF > 1/3 & BF <= 1) & (d > .2 & d <= .5))
  d.b5.2 = subset(d.bf, (BF > 1 & BF <= 3) & (d > .2 & d <= .5))
  d.b4.2 = subset(d.bf, (BF > 3 & BF <= 10) & (d > .2 & d <= .5))
  d.b3.2 = subset(d.bf, (BF > 10 & BF <= 30) & (d > .2 & d <= .5))
  d.b2.2 = subset(d.bf, (BF > 30 & BF <= 100) & (d > .2 & d <= .5))
  d.b1.2 = subset(d.bf, (BF > 100) & (d > .2 & d <= .5))
  
  
  d.b7.1 = subset(d.bf, (BF > .1 & BF <= 1/3) & (d <= .2))
  d.b6.1 = subset(d.bf, (BF > 1/3 & BF <= 1) & (d <= .2))
  d.b5.1 = subset(d.bf, (BF > 1 & BF <= 3) & (d <= .2))
  d.b4.1 = subset(d.bf, (BF > 3 & BF <= 10) & (d <= .2))
  d.b3.1 = subset(d.bf, (BF > 10 & BF <= 30) & (d <= .2))
  d.b2.1 = subset(d.bf, (BF > 30 & BF <= 100) & (d <= .2))
  d.b1.1 = subset(d.bf, (BF > 100) & (d <= .2))
  
    
  db.cat.sizes = lapply(mget(ls(pattern = "d\\.b\\d")), nrow)
  total.db.cat.sizes = Data.size
  
  
  bf.Decisive = nrow(subset(d.bf, BF > 100))
  bf.very.strong = nrow(subset(d.bf, (BF > 30) & (BF <= 100)))
  bf.strong = nrow(subset(d.bf, (BF > 10) & (BF <= 30)))
  bf.substan.H1 = nrow(subset(d.bf, (BF > 3) & (BF <= 10)))
  bf.Anecd.H1 = nrow(subset(d.bf, (BF > 1) & (BF <= 3)))
  bf.Anecd.H0 = nrow(subset(d.bf, (BF > 1/3) & (BF <= 1)))
  bf.substan.H0 = nrow(subset(d.bf, (BF > .1) & (BF <= 1/3)))
  
  bf.cat.check = sum(bf.Decisive, bf.very.strong, bf.strong, bf.substan.H1, bf.Anecd.H1, bf.Anecd.H0, bf.substan.H0)
  
  
  d.small = nrow(subset(d.bf, d <= .2)) 
  d.smtomed = nrow(subset(d.bf, (d > .2 & d <= .5)))
  d.medtolg = nrow(subset(d.bf, (d > .5 & d <= .8)))
  d.large = nrow(subset(d.bf, d > .8))
  
  d.cat.check = sum(d.small, d.smtomed, d.medtolg, d.large)
  
  
  result = matrix(unlist(db.cat.sizes), ncol = 4, byrow = TRUE)
  dimnames(result) = list(BF = c("Decisive", "Very Strong", "Strong", "Substantial", "Anecdotal", "Anecdotal(H0)", "Substantial(H0)")
                          ,"Cohen's d" = c("Small", "Small-med", "Med-large", "Large"))
  
  d.marginal = matrix(c(d.small, d.smtomed, d.medtolg, d.large), nrow = 1)
  bf.marginal = matrix(c(bf.Decisive, bf.very.strong, bf.strong, bf.substan.H1,
                         bf.Anecd.H1, bf.Anecd.H0, bf.substan.H0, bf.cat.check), ncol = 1)
  
  X.sq = chisq.test(result, correct = FALSE)
  
  chisq.table = cbind(rbind(result, d.marginal), bf.marginal)
  dimnames(chisq.table) = list(BF = c("Decisive", "Very Strong", "Strong", "Substantial", "Anecdotal", "Anecdotal(H0)", "Substantial(H0)", "Marginal.d")
                               ,"Cohen's d" = c("Small", "Small-med", "Med-large", "Large", "Marginal.BF"))
  
  
  
  # Regular expression to create all 28 plot names:
  plot.names = noquote(sprintf("d.b%d.%d", 1:7, rep(1:4, each = 7)))
  
  
  # Split screen:
  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))
  
  png("Plot_d_BF.png", res = 500, width = 4.5, height = 5.3, units = "in")
  par(mai = rep(.5, 4)) 
  par(mfcol = c(7, 4), mar = rep(.07, 4), oma = c(4.5, 5, 4.2, 5.1)) 
  
 # Only for graphical device: par(mfcol = c(7, 4), mar = rep(.08, 4), oma = rep(7, 4))
  
  
  ## Multiple plotting:
  invisible(lapply(mget(plot.names), function(db) 
    if(nrow(db) == 0){
      plot.new(); box()
    }else{
      plot(db, pch = 21, bg = 3, cex = 1.3, xaxt = "n", yaxt = "n")
    }))
  

  
  ## Advanced Labeling:
  
  mid.x = seq(grconvertX(0 + (1 /  8), "nic"), grconvertX(1 - (1 /  8), "nic"), l = 4)
  mid.y = seq(grconvertY(0 + (1 / 14), "nic"), grconvertY(1 - (1 / 14), "nic"), l = 7)
  
  gap.x = sort(c(mid.x[-length(mid.x)] + diff(mid.x)[1L] / 2, grconvertX(0:1, "nic")))
  gap.y = sort(c(mid.y[-length(mid.y)] + diff(mid.y)[1L] / 2, grconvertY(0:1, "nic")))
  
  xlim = c( grconvertX(0, "nic"), grconvertX(1, "nic") )
  ylim = c( grconvertY(0, "nic"), grconvertY(1, "nic") )
  
  
  v = c(Substantial = bf.substan.H0, Anecdotal = bf.Anecd.H0, Anecdotal = bf.Anecd.H1, Substantial = bf.substan.H1, 
        Strong = bf.strong, "Very Strong" = bf.very.strong, Decisive = bf.Decisive)/Data.size*1e2
  l = paste0(names(v), "\n", v, "%")       
  
  mtext(l, side = 4, at = mid.y, las = 1, line = .2, cex = .7, font = 2)
  

  v = paste0(c(d.small, d.smtomed, d.medtolg, d.large)/Data.size*1e2, "%")
  l = c("Small", "Medium", "Large")
  
  text(gap.x[2:4], ylim[2] + .35, l, xpd = NA, cex = 1.1, font = 2)
  text(mid.x, ylim[2] + .1, v, xpd = NA, cex = 1.1, font = 2)
  
  
  l = c('0', '.2', '.5', '.8', '> 5') 
  mtext(l, side = 1, at = gap.x, line = .55, cex = .9, font = 2)
  
  
  l = c('1/10', '1/3', '1', '3', '10', '30', '100', '> 100000')
  text(xlim[1], gap.y, l, xpd = NA, adj = 1.05, cex = 1.3, font = 2)
  
  mtext("Effect Size", side = 1, at = mean(gap.x), line = 3.5, cex = 1.2, font = 2)
  text(xlim[1], mean(gap.y), "Bayes Factor", xpd = NA, adj = c(.5, -2.2), cex = 1.8, font = 2, srt = 90)
  
  dev.off()

  
  list(chisq.table = chisq.table, X.sq = X.sq)
  
}




p.value_d = function(){
  
  ## Subsetting for a 4x4 plot p.value against BF:
  p.d4.4 = subset(pval.d, (d <= .2) & (p > .05) )
  p.d3.4 = subset(pval.d, (d > .2 & d <= .5) & (p > .05) )
  p.d2.4 = subset(pval.d, (d > .5 & d <= .8) & (p > .05) )
  p.d1.4 = subset(pval.d, (d > .8) & (p > .05) )
  
  p.d4.3 = subset(pval.d, (d <= .2) & (p > .01 & p <= .05) )
  p.d3.3 = subset(pval.d, (d > .2 & d <= .5) & (p > .01 & p <= .05) )
  p.d2.3 = subset(pval.d, (d > .5 & d <= .8) & (p > .01 & p <= .05) )
  p.d1.3 = subset(pval.d, (d > .8) & (p > .01 & p <= .05) )
  
  p.d4.2 = subset(pval.d, (d <= .2) & (p > .001 & p <= .01) )
  p.d3.2 = subset(pval.d, (d > .2 & d <= .5) & (p > .001 & p <= .01) )
  p.d2.2 = subset(pval.d, (d > .5 & d <= .8) & (p > .001 & p <= .01) )
  p.d1.2 = subset(pval.d, (d > .8) & (p > .001 & p <= .01) )
  
  p.d4.1 = subset(pval.d, (d <= .2) & (p <= .001) )
  p.d3.1 = subset(pval.d, (d > .2 & d <= .5) & (p <= .001) )
  p.d2.1 = subset(pval.d, (d > .5 & d <= .8) & (p <= .001) )
  p.d1.1 = subset(pval.d, (d > .8) & (p <= .001) )
  

  
  pd.cat.sizes = lapply(mget(ls(pattern = "p\\.d\\d")), nrow)
  total.pd.cat.sizes = Data.size
  
  p.Decisive = nrow(subset(pval.d, p <= .001)) 
  p.substan = nrow(subset(pval.d, (p > .001 & p <= .01)))
  p.positive = nrow(subset(pval.d, (p > .01 & p <= .05)))
  p.no = nrow(subset(pval.d, p > .05))
  
  p.cat.check = sum(p.Decisive, p.substan, p.positive, p.no)
  
  
  d.small = nrow(subset(pval.d, d <= .2)) 
  d.smtomed = nrow(subset(pval.d, (d > .2 & d <= .5)))
  d.medtolg = nrow(subset(pval.d, (d > .5 & d <= .8)))
  d.large = nrow(subset(pval.d, d > .8))
  
  d.cat.check = sum(d.small, d.smtomed, d.medtolg, d.large)
  

  result = matrix(unlist(pd.cat.sizes), ncol = 4, byrow = TRUE)
  rownames(result) = c("Large", "Med-large", "Small-med", "Small")
  colnames(result) = c("Decisive", "Substantial", "Positive", "None")

  d.marginal = matrix(c(d.large, d.medtolg, d.smtomed, d.small), ncol = 1)
  p.marginal = matrix(c(p.Decisive, p.substan, p.positive, p.no, p.cat.check), nrow = 1)
  
  
  X.sq = chisq.test(result, correct = FALSE)
  
  chisq.table = rbind(cbind(result, d.marginal), p.marginal)
  dimnames(chisq.table) = list("Cohen's d" = rev(c("Marginal.p", "Small", "Small-med", "Med-large", "Large")),
                               p.value = c("Decisive", "Substantial", "Positive", "None", "Marginal.d"))
  
  
  plot.names = noquote(sprintf("p.d%d.%d", 1:4, rep(1:4, each = 4)))
  

  # Split screen:
  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))
  
  png("Plot_P_d.png", res = 500, width = 4.5, height = 5.3, units = "in")
  par(mai = rep(.5, 4)) 
  par(mfcol = c(4, 4), mar = rep(.07, 4), oma = c(4.5, 5, 4.2, 5.1)) 
  
  # Only for graphical device: par(mfcol = c(4, 4), mar = rep(.08, 4), oma = rep(7, 4))
  
  
  ## Multiple plotting:
  invisible(lapply(mget(plot.names), function(pd) 
    if(nrow(pd) == 0){
      plot.new(); box()
    }else{
      plot(pd, pch = 21, bg = 3, cex = 1.3, xaxt = "n", yaxt = "n")
    }))
  
  
  
  ## Advanced Labeling:
  
  mid.x = seq(grconvertX(0 + (1 /  8), "nic"), grconvertX(1 - (1 /  8), "nic"), l = 4)
  mid.y = seq(grconvertY(0 + (1 /  8), "nic"), grconvertY(1 - (1 /  8), "nic"), l = 4)
  
  gap.x = sort(c(mid.x[-length(mid.x)] + diff(mid.x)[1L] / 2, grconvertX(0:1, "nic")))
  gap.y = sort(c(mid.y[-length(mid.y)] + diff(mid.y)[1L] / 2, grconvertY(0:1, "nic")))
  
  xlim = c( grconvertX(0, "nic"), grconvertX(1, "nic") )
  ylim = c( grconvertY(0, "nic"), grconvertY(1, "nic") )
  
  
  v = paste0(c(d.small, d.smtomed, d.medtolg, d.large)/Data.size*1e2,"%")
  l = c("Large", "Medium", "Small")       
  
  mtext(v, side = 4, at = mid.y, las = 1, line = .2, xpd = NA, cex = .75, font = 2)
  mtext(l, side = 4, at = gap.y[4:2], las = 1, line = .2, xpd = NA, cex = .8, font = 2)

  
  v = c(Decisive= p.Decisive, Substantial= p.substan, Positive= p.positive, None= p.no)/Data.size*1e2
  l = paste0(names(v), "\n", v, "%")
  
  text(mid.x, ylim[2] + .027, l, xpd = NA, cex = 1.1, font = 2)
  
  
  l = c('0', '.001', '.01', '.05', '1') 
  mtext(l, side = 1, at = gap.x, line = .55, cex = .9, font = 2)
  
  
  l = c('0', '.2', '.5', '.8', '> 5') 
  text(xlim[1], gap.y, l, xpd = NA, adj = 1.2, cex = 1.4, font = 2)
  
  
    mtext(expression(bolditalic("p")*bold("-value")), side = 1, at = mean(gap.x), line = 3.5, cex = 1.2)
    text(xlim[1], mean(gap.y), "Effect Size", xpd = NA, adj = c(.5, -2.2), cex = 1.8, font = 2, srt = 90)
    text(mean(gap.x), ylim[2] + .09, bquote(bold("Evidence against"~ bolditalic(H)[0])), cex = 1.8, xpd = NA)
    
    
  dev.off()
  
  list(chisq.table = chisq.table, X.sq = X.sq)
  
}
