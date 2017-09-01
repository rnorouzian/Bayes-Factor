type.s.mfun = function(n1, n2 = NA, d.min = 0, d.max = 1.4, alpha = .05){
  
   type.s.m = function(n1, n2, d, alpha){
    
          N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
         df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
       d.SE = 1/sqrt(N) ; ncp = d*sqrt(N)
         CI = qt(c(alpha/2, 1-(alpha/2)), df)*d.SE
    
type.s.area = pt(ifelse(d > 0, CI[1]/d.SE, CI[2]/d.SE), df, ncp, lower.tail = ifelse(d > 0, TRUE, FALSE))
      power = type.s.area + pt(ifelse(d > 0, CI[2]/d.SE, CI[1]/d.SE), df, ncp, lower.tail = ifelse(d > 0, FALSE, TRUE))
     type.s = type.s.area / power
   random.d = rt(1e4, df, ncp)*d.SE
        sig = if(d > 0) abs(random.d) > CI[2] else -abs(random.d) < CI[1]
exaggration = if(d > 0) mean(abs(random.d)[sig])/ d else mean(-abs(random.d)[sig])/ d
    
    list(exaggration = exaggration, type.s = type.s, power = power)
  }
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), las = 1)  
  
       d_range = seq(d.min, d.max, by = 5e-3)
             n = length(d_range)
         power = numeric(n)
        type.s = numeric(n)
   exaggration = numeric(n)
  
  for (i in 1L:n){
             a = type.s.m(d = d_range[i], n1 = n1, n2 = n2, alpha = alpha)
      power[i] = a$power
     type.s[i] = a$type.s
exaggration[i] = a$exaggration
  }
  
  plot(power, type.s, ty = "l", xaxt = "n", lwd = 2, font.lab = 2, col = 2)
  axis(1, at = c(alpha, seq(.2, 1, by = .2)))
  abline(v = alpha, col = 8)
  plot(power, exaggration, ty = "l", ylim = c(1, 10), xaxt = "n", yaxt = "n", lwd = 2, font.lab = 2, col = 4)
  axis(1, at = c(alpha, seq(.2, 1, by = .2)))
  axis(2, at = seq(1, 10, by = 2))
  abline(h = 1, v = alpha, col = 8)
  
}
# Example of use:
type.s.mfun(n1 = 10, alpha = .05)
