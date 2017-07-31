
type.s.m = function(n1 = 20, n2 = NA, d = .1, obs.d = .6){
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), xpd = TRUE)  
  
     N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  d.SE = 1/sqrt(N) ; ncp = d*sqrt(N)
  
 min.d = qt(1e-4, df)*d.SE  ;  max.d = qt(0.9999, df)*d.SE
  
`d|H0` = curve( dt(x/d.SE, df)/d.SE, min.d, max.d, n = 1e4, xlab = "Effect Size", 
                  ylab = NA, font = 2, font.lab = 2, ty = "n", yaxt = "n", bty = "n",
                  cex.axis = 1, cex.lab = 1, yaxs = "i")
  
    CI = qt(c(.025, .975), df)*d.SE
  
     x = seq(min.d, CI[1], l = 1e4) ;  y = dt(x /d.SE, df)/d.SE
    xx = seq(max.d, CI[2], l = 1e4) ; yy = dt(xx/d.SE, df)/d.SE
  
  polygon(c(min.d,  x, CI[1]), c( y[1],  y, rev( y[1])), col = 2, border = NA)  
  polygon(c(max.d, xx, CI[2]), c(yy[1], yy, rev(yy[1])), col = 2, border = NA)  
  
  lines(`d|H0`, lwd = 2)
  
  points(obs.d, 0, pch = 23, bg = 3, cex = 1.4, xpd = TRUE)
  
  legend("topright", "Observed \nS.S. Effect", pch = 23, pt.bg = 3, pt.cex = 1.2, bty = "n", text.font = 2)   
  abline(v = 0, col = 2, xpd = FALSE) 
  
  par(mar = c(5, 4, 1, 2))
  
`d|H1` = curve( dt(x/d.SE, df, ncp)/d.SE, min.d, max.d, n = 1e4, xlab = "Effect Size", 
                  ylab = NA, font = 2, font.lab = 2, yaxt = "n", bty = "n",
                  cex.axis = 1, cex.lab = 1, yaxs = "i", ty = "n")
  
     x = seq(min.d, CI[1], l = 1e4)   ;  y = dt(x /d.SE, df, ncp)/d.SE
    xx = seq(max.d, CI[2], l = 1e4)   ; yy = dt(xx/d.SE, df, ncp)/d.SE 
  
  polygon(c(min.d,  x, CI[1]), c( y[1],  y, rev( y[1])), col = 2, border = NA)  
  polygon(c(max.d, xx, CI[2]), c(yy[1], yy, rev(yy[1])), col = 2, border = NA) 
  
  lines(`d|H1`, lwd = 2)
  
  axis(1, at = d, col = 4, col.axis = 4, font = 2)
  points(obs.d, 0, pch = 23, bg = 3, cex = 1.4, xpd = TRUE)
  abline(v = d, col = 4, xpd = FALSE)
  
  segments(c(CI[1], CI[2]), 0, c(CI[1], CI[2]), 20, lty = 2, col = 2, xpd = NA)
  
  type.s.area = pt(ifelse(d > 0, CI[1]/d.SE, CI[2]/d.SE), df, ncp, lower.tail = ifelse(d > 0, TRUE, FALSE))
        power = type.s.area + pt(CI[2]/d.SE, df, ncp, lower.tail = FALSE)
       type.s = type.s.area / power
      p.value = 2*pt(abs(obs.d)/d.SE, df, lower.tail = FALSE)
     random.d = rt(n = 1e6, df, ncp)*d.SE
          sig = abs(random.d) > CI[2]
  exaggration = mean(abs(random.d)[sig]) / d
  
  list(exaggration.ratio = exaggration, type.s = type.s, power = power, Crit.d = CI[2], p.value = p.value)
}

# Example used in the Supplementry Document:

 type.s.m(n1 = 20, d = .1, obs.d = .6) 
