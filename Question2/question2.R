#initial starting concentrations and k constant values
y0 <- c(1, 10, 0, 0, 0) 
start_c <- c(100, 600, 150)

#function that contains equation to be solved, 4x equation and parameters - k1,k2,k3
f <- function(t, y, start_c) {
  
  #k constants
  k1 <- start_c[1] #uM/min
  k2 <- start_c[2] #uM/min
  k3 <- start_c[3] #uM/min
  
  #concentations
  E <- y[1] #uM
  S <- y[2] #uM
  ES <- y[3] #uM
  P <- y[4] #uM
  
  #derived mass law equations
  dEdt <- -k1*E*S + k2*ES + k3*ES
  dSdt <- -k1*E*S + k2*ES
  dESdt <- k1*E*S - k2*ES - k3*ES
  dPdt <- k3*ES
  
  #expressing v as a function of dP/dS
  #dPdt <- k3*ES, sub/w dESdt <- k1*E*S - k2*ES - k3*ES
  v <- k3*(k1*E*S - k2*ES - k3*ES)
  return(list(c(dEdt, dSdt, dESdt, dPdt, v)))
  
}

#total number of steps (y-axis time) and step size (step size = dt), n_steps = total number of steps for calculation
dt <- 0.001
n_steps <- 500

#start blank array (stored as y)
y <- matrix(0, nrow = n_steps + 1, ncol = 5)
t <- numeric(n_steps + 1)  #t in minutes
y[1,] <- y0 #start position to initial starting concentration defined above
t[1] <- 0 #t in minutes, start at 0
colnames(y) <- c("E", "S", "ES", "P", "v")

#calculation with steps in loop
for (i in 1:n_steps) {
  #find k1
  k1 <- f(t[i], y[i,], start_c) #function(x0+1, y0+1, start conc.)
  #find k2 using k1
  k2 <- f(t[i] + 0.5*dt, y[i,] + 0.5*dt*k1[[1]], start_c) #function(x0+1, y0+1, start conc.)
  #find k3 using k2
  k3 <- f(t[i] + 0.5*dt, y[i,] + 0.5*dt*k2[[1]], start_c)
  #find k4 using k3
  k4 <- f(t[i] + dt, y[i,] + dt*k3[[1]], start_c)
  #average of the gradients
  y[i+1,] <- y[i,] + dt*(k1[[1]] + 2*k2[[1]] + 2*k3[[1]] + k4[[1]])/6
  #tn <- +1 to time
  t[i+1] <- t[i] + dt
  
}

#Plot values (stored as y, against t (time))
plot(t, y[,2], type="l", xlab="Time (min)", ylab="Concentration (μM)", main="Enzyme Kinetic Model", ylim=(c(-0.5, 10)), xlim=(c(-0, 0.7)))
lines(t, y[,1], col="orange") #E
lines(t, y[,3], col="red") #ES
lines(t, y[,4], col="blue") #P
legend("right", legend=c("S", "ES", "P", "E"), col=c("black", "red", "blue", "orange"), lty=1)

#plot v as f[S]
plot(y[,2], y[,5], col="purple", type="l", xlab="Substrate (μM)", ylab="rate (μM/min)", main="Product as f(Substrate)", ylim=(c(0, 100)), xlim=(c(-0, 10)))
legend("right", legend=c("v"), col=c("purple"), lty=1)
Vmax <- max(y[,5])
Vmax_Sindex <- match(Vmax, y[,5])
Vmax_S <- y[Vmax_Sindex,2]
paste("Vmax is", Vmax, "uM/min when [S] is at", Vmax_S, "uM")
