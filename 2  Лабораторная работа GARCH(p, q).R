# 2 Лабораторная работа: GARCH(p, q)
library("tseries") 
# Построить график стационарного процесса {ℎ𝑛} и график волатильности {𝜎𝑛} процесса 𝐺𝐴𝑅𝐶𝐻(1,0), из 𝑛 = 1000 наблюдений.
n = 1000
# Стационарность достигается при 0 < 𝑎i < 1
a_0 = 0.7; a_1 = 0.3 
h_0 = 0.5
b_0 = 0.2; b_1 = 0.5

garch = function(n,h0,a0,a1,ret) 
{h = array(dim = n); x = array(dim = n)
 E = rnorm(n)
 x[1] = a0 + a1*h0^2
 h[1] = sqrt(x[1])*E[1]
 for (i in 2:n) {x[i] = a0 + a1*h[i-1]^2
                 h[i] = sqrt(x[i])*E[i]}
 if (ret == 1) return(h) else return(sqrt(x))
}

par(mfrow = c(2, 1))
garch1 = garch(n,h_0,a_0,a_1,1)
volat = garch(n,h_0,a_0,a_1,0)

plot(garch1, type = 'l', col = 'purple', main = "график стационарного процесса GARCH(1,0)")
plot(volat, type = 'l', col = 'red', main = "график волатильности процесса GARCH(1,0)")

#Оценка параметров a0 и a1 процесса AR(1) с помощью МНК
MNK_0 = function(n,h,a_1) {sum = 0
                           for (i in 2:n)
                           sum = sum + (h[i]^2 - a_1*h[i-1]^2)
                           return(sum/n)}
MNK_1 = function(n,h,a_0) {sum_1 = 0; sum_2 = 0
                           for (i in 2:n) {sum_1 = sum_1 + h[i-1]^2*(h[i]^2 - a_0)
                                           sum_2 = sum_2 + h[i-1]^4}
                           return (sum_1/sum_2)}
MNK_0 = MNK_0(n,garch1,a_1);MNK_0
MNK_1 = MNK_1(n,garch1,a_0);MNK_1

#Оценка параметров a0 и a1 процесса AR(1) при помощи функции garch():
tseries::garch(garch1,order=c(1,0),start=c(a_0,a_1))

#Построить стационарный процесс 𝐺𝐴𝑅𝐶𝐻(1,1), из 𝑛 = 1000 наблюдений и оценить его параметры
garch_11 = function(n,h0,b0,b1,a0,a1,ret) 
{h = array(dim = n); x = array(dim = n)
 E = rnorm(n)
 x[1] = a0 + a1*(h0)^2+b1*(b0)^2
 h[1] = sqrt(x[1])*E[1]
for (i in 2:n) {x[i] = a0 + a1*h[i-1]^2 + b1*x[i-1]
                h[i] = sqrt(x[i])*E[i]}
if (ret == 1) return(h) else return(sqrt(x))
}
garch11 = garch_11(n,h_0,b_0,b_1,a_0,a_1,1)
tseries::garch(garch11,order=c(1,1),start=c(a_0,a_1,b_1))
plot(garch11, type = 'l', col = 'purple', main = "график стационарного процесса GARCH(1,1)")

#Построить график стационарного процесса 𝐺𝐴𝑅𝐶𝐻(3,0), из 𝑛 = 1100 наблюдений.
n = 1100
a = c(0.9,0.8,0.3,0.1)
h_0 = c(1,0.5,0.9)
n_1 = 1000; n_2 = 100

garch3 = function(n,h0,a,ret) 
{h = array(dim = n); x = array(dim = n)
 E = rnorm(n,0,1)
 x[1] = a[1] + a[2]*h0[3]^2 + a[3]*h0[2]^2 + a[4]*h0[1]^2
 h[1] = sqrt(x[1])*E[1]
 x[2] = a[1] + a[2]* h[1]^2 + a[3]*h0[3]^2 + a[4]*h0[2]^2
 h[2] = sqrt(x[2])*E[2]
 x[3] = a[1] + a[2]* h[2]^2 + a[3]* h[1]^2 + a[4]*h0[3]^2
 h[3] = sqrt(x[3])*E[3]
  
for (i in 4:n) {x[i] = a[1] + a[2]*h[i-1]^2 + a[3]*h[i-2]^2 + a[4]*h[i-3]^2
                h[i] = sqrt(x[i])*E[i]}
if (ret == 1) return(h) 
else return(sqrt(x))
 }

garch_3 = garch3(n,h_0,a,1)
h_learn = array(dim = n_1)

for(i in 1:n_1){h_learn[i] = garch_3[i]}

tseries::garch(h_learn,order=c(3,0),start=a)$coef
ap = tseries::garch(h_learn,order=c(3,0),start=a)$coef

p = function(n_1,n_2,h_0,h,a) 
{hp = array(dim=n_1+n_2)
   for(i in 1:n_1) {hp[i] = h[i]^2
                    hp[n_1+1] = a[1] + a[2]*hp[n_1]   + a[3]*h[n_1-1]^2 + a[4]*h[n_1-2]^2
                    hp[n_1+2] = a[1] + a[2]*hp[n_1+1] + a[3]*h[n_1]^2   + a[4]*h[n_1-1]^2
                    hp[n_1+3] = a[1] + a[2]*hp[n_1+2] + a[3]*hp[n_1+1]  + a[4]*h[n_1]^2}
  
  for(i in (n_1+4):(n_1+n_2)) {hp[i] = a[1] + a[2]*hp[i-1] + a[3]*hp[i-2] + a[4]*hp[i-3]
                               hp = sqrt(hp)}
  for(i in 1:n_1){hp[i] = h[i]}
return(hp)
}

hp2 = p(n_1,n_2,h0,h_learn,ap)

par(mfrow = c(2, 1))
plot(garch_3, type='l',col='purple', main = "график стационарного процесса GARCH(3,0)")
plot(garch_3, type='l',col='red', main = "график  процесса GARCH(3,0)+прогнозы")
lines(hp2, type='l', col='blue')
