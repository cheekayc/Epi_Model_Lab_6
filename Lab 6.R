library(tidyverse)
library(deSolve)

# 1. Setting up the Ebola model

EBOLA = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
  
    # rate of change
    dS =  -S/N * (beta.I * I + beta.H * H + beta.F * FF);
    dE = S/N * (beta.I * I + beta.H * H  +  beta.F * FF)  -  alpha * E
    dI = alpha * E - I * (gamma.h * theta1 + gamma.d * (1 - theta1) * delta1 + gamma.i * (1 - theta1) * (1 - delta1));
    dH = gamma.h * theta1 * I - gamma.dh * delta2 * H  -  gamma.ih * (1 - delta2) * H
    dFF =  gamma.dh * delta2 * H  +  gamma.d * (1 - theta1) * delta1 * I  -  gamma.f * FF
    dR =  gamma.f * FF  +  gamma.i * (1 - theta1) * (1 - delta1) * I  +  gamma.ih * (1 - delta2) * H
    
    dcumInci = alpha * E # cumulative incidence
    
    # return the rate of change
    list(c(dS,  dE,  dI,  dH,  dFF,  dR,  dcumInci))
  }) # end with(as.list...)
}

####################################################################################
## Simulate the DRC 1995 outbreak
####################################################################################

# Initial conditions:
N = 2e5; 
E0 = H0 = FF0 = R0 = 0; 
I0 = 3; 
cumInci0 = 3; 
S0 = N - I0;
state = c(S = S0, E = E0, I = I0, H = H0, FF = FF0, R = R0, cumInci = cumInci0);

## parameters for the DRC 1995 outbreak (Tables 3 & 4 in Legrand et al. 2007)
# NO INTERVENTION
alpha = 1/7; # incubation period: 7 days
gamma.h = 1/5; # from onset to hospitalization: 5 days
gamma.d = 1/9.6; # from onset to death: 9.6 days
gamma.i = 1/10; # from onset to end of infectiousness for survivors: 10 days
gamma.f = 1/2; # from death to traditional burial 2 days
gamma.ih = 1/(10 - 5); # from hospitalization to end of infectiousness for survivors
gamma.dh = 1/(9.6 - 5); # from hospitalization to death
theta1 = 0.67; # proportion infectious in the hospital
delta1 = 0.8; # CFR for non-hospitalized
delta2 = 0.8; # CFR for hospitalize
p = 0.85; ## reduced by a factor of p as we are running it deterministically
beta.I = 0.588/7; # transmission rate in the community
beta.H = 0.794/7; # transmission rate in the hospital
beta.F = 7.653/7 * p # transmission rate at funerals; 

# Put all the parameters together: parms without intervention
parmsNoCtrl = c(alpha = alpha, # incubation period: 7 days
              gamma.h = gamma.h, # from onset to hospitalization: 5 days
              gamma.d = gamma.d, # from onset to death: 9.6 days
              gamma.i = gamma.i, # from onset to end of infectiousness for survivors: 10 days
              gamma.f = gamma.f, # from death to traditional burial 2 days
              gamma.ih = gamma.ih, # from hospitalization to end of infectiousness for survivors
              gamma.dh = gamma.dh, # from hospitalization to death
              theta1 = theta1, # proportion infectious in the hospital
              delta1 = delta1, # CFR for non-hospitalized
              delta2 = delta2, # CFR for hospitalize
              beta.I = beta.I, # transmission rate in the community
              beta.H = beta.H, # transmission rate in the hospital
              beta.F = beta.F # transmission rate at funerals, 
);

# simulation time steps
times = seq(1, 140, by = 1)

# Call the function
simNoCtrl = ode(y = state, times = times, func = EBOLA, parms = parmsNoCtrl)

# Get weekly incidence
cuminci_week = data.frame(matrix(ncol = 3, nrow = 20))
colnames(cuminci_week) = c("week", "cum_incidence", "weekly_incidence")

for (i in 1:20) {
  cuminci_week[i, 'week'] = i
  cuminci_week[i, 'cum_incidence'] = simNoCtrl[simNoCtrl[ , 'time'] == 7*i, 'cumInci']
  if(i == 1) {cuminci_week[i, 'weekly_incidence'] = cuminci_week[i, 'cum_incidence']} 
  else {cuminci_week[i, 'weekly_incidence'] = cuminci_week[i, 'cum_incidence'] - cuminci_week[i-1, 'cum_incidence']}
}

# If no intervention were in place, what would happen? Plot the results.
plot(x = cuminci_week[, 1], y = cuminci_week[, 3], xlab = 'Time (week)', ylab = 'New incidence per week', main = 'No intervention: Weekly Incidence vs. Time', type = 'l')

#######################################
## NOW RUN IT WITH INTERVENTION AFTER 9 WEEK OF FIRST CASE
#######################################

# Parameters with intervention
z = 0.88; # effectiveness of intervention in the community
z.H = 0.75; # effectiveness of intervention in the hospital
z.F = 0.875; # effectiveness of intervention in safe burial
parmsCtrl = c(alpha = alpha, # incubation period: 7 days
              gamma.h = gamma.h, # from onset to hospitalization: 5 days
              gamma.d = gamma.d, # from onset to death: 9.6 days
              gamma.i = gamma.i, # from onset to end of infectiousness for survivors: 10 days
              gamma.f = gamma.f, # from death to traditional burial 2 days
              gamma.ih = gamma.ih, # from hospitalization to end of infectiousness for survivors
              gamma.dh = gamma.dh, # from hospitalization to death
              theta1 = theta1, # proportion infectious in the hospital
              delta1 = delta1, # CFR for non-hospitalized
              delta2 = delta2, # CFR for hospitalize
              beta.I = beta.I * (1-z), # transmission rate in the community after intervention
              beta.H = beta.H * (1-z.H), # transmission rate in the hospital after intervention
              beta.F = beta.F * (1-z.F) # transmission rate at funerals after intervention 
);

## simulation starts on Mar 1, no intervention until May 4 (~9 weeks later)
## last case identified on July 12 (~19w weeks later)
times1 = seq(1, 9*7, by = 1); # first 9 weeks with out intervention
times2 = seq(9*7, 20*7, by = 1); # afterward with incidence

# Stage 1: no intervention
simNoCtrl = ode(y = state, times = times1, func = EBOLA, parms = parmsNoCtrl);

# Stage 2: yes intervention
## Use the numbers in the last day in Stage 1 (the 63rd day) as the initial conditions for Stage 2:
state2 = c(S = tail(simNoCtrl[ , 'S'], 1), E = tail(simNoCtrl[ , 'E'], 1),
         I = tail(simNoCtrl[ , 'I'], 1), H = tail(simNoCtrl[ , 'H'], 1),
         FF = tail(simNoCtrl[ , 'FF'], 1), R = tail(simNoCtrl[ , 'R'], 1),
         cumInci = tail(simNoCtrl[ , 'cumInci'] ,1));

simCtrl = ode(y = state2, times = times2, func = EBOLA, parms = parmsCtrl);

# Combine the two stages. Remove the redundant 63rd day for Stage 2 during combine.
sim2 = rbind(simNoCtrl, simCtrl[-1, ]) # the [-1, ] for the simCtrl is telling R to delete the the 1st row of simCtrl during rbind

# Get weekly incidence
week_inci_2 = data.frame(matrix(ncol = 3, nrow = 20))
colnames(week_inci_2) = c("week", "cum_incidence", "weekly_incidence")

for (i in 1:20) {
  week_inci_2[i, 'week'] = i
  week_inci_2[i, 'cum_incidence'] = sim2[sim2[ , 'time'] == 7*i, 'cumInci']
  if(i == 1) {week_inci_2[i, 'weekly_incidence'] = week_inci_2[i, 'cum_incidence']} 
  else {week_inci_2[i, 'weekly_incidence'] = week_inci_2[i, 'cum_incidence'] - week_inci_2[i-1, 'cum_incidence']}
}

# What happen when intervention measures were implemented from Day 64 onwards? Plot the results.
plot(x = week_inci_2[, 1], y = week_inci_2[, 3], xlab = 'Time (week)', ylab = 'New incidence per week', main = 'Intervention (Day 64): Weekly Incidence vs. Time', type = 'l')

#######################################
## COMPARE: NO-INTERVENTION, YES-INTERVENTION, OBSERVATIONS
#######################################
## mock data
dates = seq(as.Date('1995/3/4'), as.Date('1995/7/15'), by = 7);
obs = c(4, 2, 5, 1, 7, 6, 18, 25, 58, 40, 52, 28, 16, 19, 4, 5, 0, 0, 1, 0); # 291 cases
da = cbind(dates,obs)

## plot results：simulated (with intervention) v. observations
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = c(1.6, 0.5, 0), cex = 1.1)
## Note: In barplot, bars are not drawn at intervals 1:10, etc., but on something else. 
## The actual x-coordinates can be seen if you save your barplot object. 
df.bar = barplot(obs, col = 'white', xlab = 'Time from onset (week)', ylab = 'Weekly incidence', ylim = c(0, 80), main = 'Compare Observed with Simulation')
axis(1, at = df.bar, lab = 1:20)
lines(x = df.bar, y = inci2, lwd = 2, col = 'blue') # use 'x = df.bar' to match with weeks
lines(x = df.bar, y = cuminci_week[, 3], lwd = 2, col = 'red', lty = "dashed")
legend('topleft', c('Observed', 'Simulated: with intervention', 'Simulated: without intervention'),
       lty = c(NA, 1, 2), pch = c(0, NA, NA), col = c('black', 'blue', 'red'), cex = 0.8, bty = 'n')


############################################################################## 
## SENSITIVITY TESTS & EFFECTIVENESS OF INTERVENTION
##############################################################################
# Change the strength of intervention on transmission within the community (z) to 0.2, 0.4, 0.6, and 0.8 (keep other parms/settings the same).
## Reduction in transmission within the community:
times1 = seq(1, 9*7, by = 1); # first wk weeks with out intervention
times2 = seq(9*7, 20*7, by = 1);

zs = seq(0.2, 0.8, by = 0.2)
res = matrix(0, 20, length(zs))

for (i in 1:length(zs)){
  z = zs[i]; # effectiveness of intervention in the community
  z.H = 0.75; # effectiveness of intervention in the hospital
  z.F = 0.875; # effectiveness of intervention in safe burial
  parmsCtrl = c(alpha = alpha, # incubation period: 7 days
              gamma.h = gamma.h, # from onset to hospitalization: 5 days
              gamma.d = gamma.d, # from onset to death: 9.6 days
              gamma.i = gamma.i, # from onset to end of infectiousness for survivors: 10 days
              gamma.f = gamma.f, # from death to traditional burial 2 days
              gamma.ih = gamma.ih, # from hospitalization to end of infectiousness for survivors
              gamma.dh = gamma.dh, # from hospitalization to death
              theta1 = theta1, # proportion infectious in the hospital
              delta1 = delta1, # CFR for non-hospitalized
              delta2 = delta2, # CFR for hospitalized
              beta.I = beta.I * (1-z), # transmission rate in the community
              beta.H = beta.H * (1-z.H), # transmission rate in the hospital
              beta.F = beta.F * (1-z.F) # transmission rate at funerals, 
  );
  
simNoCtrl = ode(y = state, times = times1, func = EBOLA, parms = parmsNoCtrl);
  state2 = c(S = tail(simNoCtrl[ , 'S'], 1), E = tail(simNoCtrl[ , 'E'], 1), I = tail(simNoCtrl[ , 'I'], 1),
           H = tail(simNoCtrl[ , 'H'], 1), FF = tail(simNoCtrl[ , 'FF'], 1), R = tail(simNoCtrl[ , 'R'], 1),
           Inci = tail(simNoCtrl[ , 'cumInci'], 1));

simCtrl = ode(y = state2, times = times2, func = EBOLA, parms = parmsCtrl);
  
sim = rbind(simNoCtrl, simCtrl[-1, ])
  
# save the result before exiting the loop
res[ , i] = sim[seq(7, nrow(sim), by = 7), 'cumInci'] - c(0, sim[seq(7, nrow(sim)-7, by = 7), 'cumInci']) # get weekly incidence
} 

var.name1 = c('z = 0.2', 'z = 0.4', 'z = 0.6', 'z = 0.8')
colnames(res) = var.name1

ymax = max(res)*1.2
par(mar = c(3, 3, 1, 1), mgp = c(2, 0.5, 0), cex = 1.1)
df.bar = barplot(obs, col = 'white', xlab = 'Date of onset (week)', ylab = 'Weekly incidence', ylim = c(0, ymax), main = 'Different strength of intervention within community (z)')
axis(1, at = df.bar, lab = 1:20)
matlines(x = df.bar, y = res, lty = 1, lwd = 2, col = rainbow(ncol(res))) # use 'x = df.bar' to match with weeks
legend(0, ymax/1.2, c(paste('z = ', zs, sep = '')), lty = 1, lwd = 2, col = rainbow(ncol(res)), bty = 'n', cex = 0.8)


# Set z back to 0.88. Change the strength of intervention on transmission within the hospital (z.H) to 0.2, 0.4, 0.6, and 0.8 (keep other parms/settings the same).
## Reduction in transmission within the hospital:
times1 = seq(1, 9*7, by = 1); # first wk weeks with out intervention
times2 = seq(9*7, 20*7, by = 1);

zh = seq(0.2, 0.8, by = 0.2)
res = matrix(0, 20, length(zh))

for (i in 1:length(zh)){
  z = 0.88; # effectiveness of intervention in the community
  z.H = zh[i]; # effectiveness of intervention in the hospital
  z.F = 0.875; # effectiveness of intervention in safe burial
  parmsCtrl = c(alpha = alpha, # incubation period: 7 days
                gamma.h = gamma.h, # from onset to hospitalization: 5 days
                gamma.d = gamma.d, # from onset to death: 9.6 days
                gamma.i = gamma.i, # from onset to end of infectiousness for survivors: 10 days
                gamma.f = gamma.f, # from death to traditional burial 2 days
                gamma.ih = gamma.ih, # from hospitalization to end of infectiousness for survivors
                gamma.dh = gamma.dh, # from hospitalization to death
                theta1 = theta1, # proportion infectious in the hospital
                delta1 = delta1, # CFR for non-hospitalized
                delta2 = delta2, # CFR for hospitalized
                beta.I = beta.I * (1-z), # transmission rate in the community
                beta.H = beta.H * (1-z.H), # transmission rate in the hospital
                beta.F = beta.F * (1-z.F) # transmission rate at funerals, 
  );
  
  simNoCtrl = ode(y = state, times = times1, func = EBOLA, parms = parmsNoCtrl);
  state2 = c(S = tail(simNoCtrl[ , 'S'], 1), E = tail(simNoCtrl[ , 'E'], 1), I = tail(simNoCtrl[ , 'I'], 1),
             H = tail(simNoCtrl[ , 'H'], 1), FF = tail(simNoCtrl[ , 'FF'], 1), R = tail(simNoCtrl[ , 'R'], 1),
             Inci = tail(simNoCtrl[ , 'cumInci'], 1));
  
  simCtrl = ode(y = state2, times = times2, func = EBOLA, parms = parmsCtrl);
  
  sim = rbind(simNoCtrl, simCtrl[-1, ])
  
  # save the result before exiting the loop
  res[ , i] = sim[seq(7, nrow(sim), by = 7), 'cumInci'] - c(0, sim[seq(7, nrow(sim)-7, by = 7), 'cumInci']) # get weekly incidence
} 

var.name1 = c('0.2', '0.4', '0.6', '0.8')
colnames(res) = var.name1

ymax = max(res)*1.2
par(mar = c(3, 3, 1, 1), mgp = c(2, 0.5, 0), cex = 1.1)
df.bar = barplot(obs, col = 'white', xlab = 'Date of onset (week)', ylab = 'Weekly incidence', ylim = c(0, ymax), main = 'Different strength of intervention within hospital (z.H)')
axis(1, at = df.bar, lab = 1:20)
matlines(x = df.bar, y = res, lty = 1, lwd = 2, col = rainbow(ncol(res))) # use 'x = df.bar' to match with weeks
legend(0, ymax/1.2, c(paste('z.H = ', zs, sep = '')), lty = 1, lwd = 2, col = rainbow(ncol(res)), bty = 'n', cex = 0.8)


# Set z = 0.88 and z.H = 0.75. Change the strength of intervention on transmission at funerals (z.F) to 0.2, 0.4, 0.6, and 0.8 (keep other parms the same).
## Reduction in transmission at funerals:
times1 = seq(1, 9*7, by = 1); # first wk weeks with out intervention
times2 = seq(9*7, 20*7, by = 1);

zf = seq(0.2, 0.8, by = 0.2)
res3 = matrix(0, 20, length(zf))

for (i in 1:length(zf)){
  z = 0.88; # effectiveness of intervention in the community
  z.H = 0.75; # effectiveness of intervention in the hospital
  z.F = zf[i]; # effectiveness of intervention in safe burial
  parmsCtrl = c(alpha = alpha, # incubation period: 7 days
                gamma.h = gamma.h, # from onset to hospitalization: 5 days
                gamma.d = gamma.d, # from onset to death: 9.6 days
                gamma.i = gamma.i, # from onset to end of infectiousness for survivors: 10 days
                gamma.f = gamma.f, # from death to traditional burial 2 days
                gamma.ih = gamma.ih, # from hospitalization to end of infectiousness for survivors
                gamma.dh = gamma.dh, # from hospitalization to death
                theta1 = theta1, # proportion infectious in the hospital
                delta1 = delta1, # CFR for non-hospitalized
                delta2 = delta2, # CFR for hospitalized
                beta.I = beta.I * (1-z), # transmission rate in the community
                beta.H = beta.H * (1-z.H), # transmission rate in the hospital
                beta.F = beta.F * (1-z.F) # transmission rate at funerals, 
  );
  
simNoCtrl = ode(y = state, times = times1, func = EBOLA, parms = parmsNoCtrl);

state2 = c(S = tail(simNoCtrl[ , 'S'], 1), E = tail(simNoCtrl[ , 'E'], 1), I = tail(simNoCtrl[ , 'I'], 1),
             H = tail(simNoCtrl[ , 'H'], 1), FF = tail(simNoCtrl[ , 'FF'], 1), R = tail(simNoCtrl[ , 'R'], 1),
             Inci = tail(simNoCtrl[ , 'cumInci'], 1));
  
simCtrl = ode(y = state2, times = times2, func = EBOLA, parms = parmsCtrl);
  
sim = rbind(simNoCtrl, simCtrl[-1, ])
  
  # save the result before exiting the loop
  res3[ , i] = sim[seq(7, nrow(sim), by = 7), 'cumInci'] - c(0, sim[seq(7, nrow(sim)-7, by = 7), 'cumInci']) # get weekly incidence
} 

var.name1 = c('0.2', '0.4', '0.6', '0.8')
colnames(res3) = var.name1

ymax = max(res3)*1.2
par(mar = c(3, 3, 1, 1), mgp = c(2, 0.5, 0), cex = 1.1)
df.bar = barplot(obs, col = 'white', xlab = 'Date of onset (week)', ylab = 'Weekly incidence', ylim = c(0, ymax), main = 'Different strength of intervention at funerals (z.F)')
axis(1, at = df.bar, lab = 1:20)
matlines(x = df.bar, y = res3, lty = 1, lwd = 2, col = rainbow(ncol(res3))) # use 'x = df.bar' to match with weeks
legend(0, ymax/1.2, c(paste('z.F = ', zf, sep = '')), lty = 1, lwd = 2, col = rainbow(ncol(res3)), bty = 'n', cex = 0.8)


# Set z = 0.88, z.H = 0.75, z.F = 0.875. Change the timing of intervention to 6, 8, 10, 12 weeks after the outbreak started (keep other parms the same)
## Simulate different intervention time:
time_var = seq(6, 12, by = 2)

res = matrix(0, 20, length(time_var))

for (i in 1:length(time_var)){
  times1 = seq(1, time_var[i]*7, by = 1); 
  times2 = seq(time_var[i]*7, 20*7, by = 1);
  z = 0.88; # intervention in the community
  z.H = 0.75; # intervention in the hospital
  z.F = 0.875; # intervention in safe burial
  parmsCtrl = c(alpha = alpha, # incubation period: 7 days
              gamma.h = gamma.h, # from onset to hospitalization: 5 days
              gamma.d = gamma.d, # from onset to death: 9.6 days
              gamma.i = gamma.i, # from onset to end of infectiousness for survivors: 10 days
              gamma.f = gamma.f, # from death to traditional burial 2 days
              gamma.ih = gamma.ih, # from hospitalization to end of infectiousness for survivors
              gamma.dh = gamma.dh, # from hospitalization to death
              theta1 = theta1, # proportion infectious in the hospital
              delta1 = delta1, # CFR for non-hospitalized
              delta2 = delta2, # CFR for hospitalized
              beta.I = beta.I * (1-z), # transmission rate in the community
              beta.H = beta.H * (1-z.H), # transmission rate in the hospital
              beta.F = beta.F * (1-z.F) # transmission rate at funerals, 
  );
  simNoCtrl = ode(y = state, times = times1, func = EBOLA, parms = parmsNoCtrl);
  state2 = c(S = tail(simNoCtrl[, 'S'], 1), E = tail(simNoCtrl[ , 'E'], 1), I = tail(simNoCtrl[ , 'I'], 1),
           H = tail(simNoCtrl[ , 'H'], 1), FF = tail(simNoCtrl[ , 'FF'], 1), R = tail(simNoCtrl[ , 'R'], 1),
           Inci = tail(simNoCtrl[ , 'cumInci'], 1));
  simCtrl = ode(y = state2, times = times2, func = EBOLA, parms = parmsCtrl);
  
  sim = rbind (simNoCtrl, simCtrl[-1, ])
  # save the result before exiting the loop
  res[ , i] = sim[seq(7, nrow(sim), by = 7), 'cumInci'] - c(0, sim[seq(7, nrow(sim) - 7, by = 7), 'cumInci']) # get weekly incidence
}

ymax = max(res) * 1.2
par(mar = c(3, 3, 1, 1), mgp = c(2, 0.5, 0), cex = 1.1)
df.bar = barplot(obs, col = 'white', xlab = 'Date of onset (week)',
               ylab = 'Weekly incidence', ylim = c(0, ymax))
axis(1, at = df.bar, lab = 1:20)
# matlines allow us to plot all the simulations at the same time
matlines(x = df.bar, y = res, lty = 1, lwd = 2, col = rainbow(ncol(res))) # use 'x = df.bar' to match with weeks
legend(0, ymax/1.2, c(paste('z = ', zs, sep = '')),
       lty = 1, lwd = 2, col = rainbow(ncol(res)), bty = 'n', cex = 0.8)

### Try on your own
times1 = seq(1, 9*7, by = 1); # first wk weeks with out intervention
times2 = seq(9*7, 20*7, by = 1);

zs = seq(0.2, 0.8, by = 0.2)
res = matrix(0, 20, length(zs))
for (i in 1:length(zs)){
  z = zs[i]; # intervention in the community
  z.H = 0.75; # intervention in the hospital
  z.F = 0.875; # intervention in safe burial
  parmsCtrl = c(alpha = alpha, # incubation period: 7 days
                gamma.h = gamma.h, # from onset to hopspitalization: 5 days
                gamma.d = gamma.d, # from onset to death: 9.6 days
                gamma.i = gamma.i, # from onset to end of infectiousness for survivors: 10 days
                gamma.f = gamma.f, # from death to traditional burial 2 days
                gamma.ih = gamma.ih, # from hospitalization to end of infectiousness for survivors
                gamma.dh = gamma.dh, # from hospitalization to death
                theta1 = theta1, # proportion infectious in the hospital
                delta1 = delta1, # CFR for unhospitalized
                delta2 = delta2, # CFR for hospitalize
                beta.I = beta.I * (1-z), # transmission rate in the community
                beta.H = beta.H * (1-z.H), # transmission rate in the hospital
                beta.F = beta.F * (1-z.F) # transmission rate at funerals, 
  );
  simNoCtrl = ode(y = state, times = times1, func = EBOLA, parms = parmsNoCtrl);
  state2 = c(S = tail(simNoCtrl[, 'S'], 1), E = tail(simNoCtrl[ , 'E'], 1), I = tail(simNoCtrl[ , 'I'], 1),
             H = tail(simNoCtrl[ , 'H'], 1), FF = tail(simNoCtrl[ , 'FF'], 1), R = tail(simNoCtrl[ , 'R'], 1),
             Inci = tail(simNoCtrl[ , 'cumInci'], 1));
  simCtrl = ode(y = state2, times = times2, func = EBOLA, parms = parmsCtrl);
  
  sim = rbind (simNoCtrl, simCtrl[-1, ])
  # save the result before exiting the loop
  res[ , i] = sim[seq(7, nrow(sim), by = 7), 'cumInci'] - c(0, sim[seq(7, nrow(sim) - 7, by = 7), 'cumInci']) # get weekly incidence
}

ymax = max(res) * 1.2
par(mar = c(3, 3, 1, 1), mgp = c(2, 0.5, 0), cex = 1.1)
df.bar = barplot(obs, col = 'white', xlab = 'Date of onset (week)',
                 ylab = 'Weekly incidence', ylim = c(0, ymax))
axis(1, at = df.bar, lab = 1:20)
# matlines allow us to plot all the simulations at the same time
matlines(x = df.bar, y = res, lty = 1, lwd = 2, col = rainbow(ncol(res))) # use 'x = df.bar' to match with weeks
legend(0, ymax/1.2, c(paste('z = ', zs, sep = '')),
       lty = 1, lwd = 2, col = rainbow(ncol(res)), bty = 'n', cex = 0.8)

