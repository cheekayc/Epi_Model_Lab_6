## Lab 6: A complex compartmental model for Ebola
## model taken from Legrand et al. 2007

#############################################################################
## !!! NOTE: Please focus on the questions in the lab slides (i.e. pdf file)  
##              -those are the core questions you need to answer. 
## Any additional questions in this R script are meant as guidelines.
#############################################################################

## if you haven't had it installed yet, type in the following first:
## install.packages('deSolve')
library(deSolve)

EBOLA=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
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
N=2e5; E0=H0=FF0=R0=0; I0=3; cumInci0=3; S0=N-I0;
state=c(S=S0,E=E0,I=I0,H=H0,FF=FF0,R=R0,cumInci=cumInci0);

## parameters for the DRC 1995 outbreak (Tables 3 & 4 in Legrand et al. 2007)
# no intervention
alpha=1/7; # incubation period: 7 days
gamma.h=1/5; # from onset to hopspitalization: 5 days
gamma.d=1/9.6; # from onset to death: 9.6 days
gamma.i=1/10; # from onset to end of infectiousness for survivors: 10 days
gamma.f=1/2; # from death to traditional burial 2 days
gamma.ih=1/(10-5); # from hospitalization to end of infectiousness for survivors
gamma.dh=1/(9.6-5); # from hospitalization to death
theta1=.67; # proportion infectious in the hospital
delta1=.8; # CFR for unhospitalized
delta2=.8; # CFR for hospitalize
p=.85; ## reduced by a factor of p as we are running it deterministically
beta.I=.588/7; # transmission rate in the community
beta.H=.794/7; # transmission rate in the hospital
beta.F=7.653/7 * p # transmission rate at funerals; 

# parms without intervention
parmsNoCtrl=c(alpha=alpha, # incubation period: 7 days
              gamma.h=gamma.h, # from onset to hopspitalization: 5 days
              gamma.d=gamma.d, # from onset to death: 9.6 days
              gamma.i=gamma.i, # from onset to end of infectiousness for survivors: 10 days
              gamma.f=gamma.f, # from death to traditional burial 2 days
              gamma.ih=gamma.ih, # from hospitalization to end of infectiousness for survivors
              gamma.dh=gamma.dh, # from hospitalization to death
              theta1=theta1, # proportion infectious in the hospital
              delta1=delta1, # CFR for unhospitalized
              delta2=delta2, # CFR for hospitalize
              beta.I=beta.I, # transmission rate in the community
              beta.H=beta.H, # transmission rate in the hospital
              beta.F=beta.F # transmission rate at funerals, 
);


####################################
## FIRST RUN WITHOUT CONTROL FOR 20 WEEKS:
####################################
# Now run the model for 20 weeks (i.e.: times=1:140) using the initial conditions on Slide #8 and the parameters for when no intervention is in place (Slides # 9)
# [Q1] If no intervention were in place, what would happen?  Plot your results (e.g. cumulative incidence or weekly incidence) and describe what you found (1 pt)

 # CODE ON YOUR OWN







#######################################
## NOW RUN IT WITH INTERVENTION AFTER 9 WEEK OF FIRST CASE
#######################################

# parms with intervention
z=.88; # intervention in the community
z.H=.75; # intervention in the hospital
z.F=0.875; # intervention in safe burial
parmsCtrl=c(alpha=alpha, # incubation period: 7 days
            gamma.h=gamma.h, # from onset to hopspitalization: 5 days
            gamma.d=gamma.d, # from onset to death: 9.6 days
            gamma.i=gamma.i, # from onset to end of infectiousness for survivors: 10 days
            gamma.f=gamma.f, # from death to traditional burial 2 days
            gamma.ih=gamma.ih, # from hospitalization to end of infectiousness for survivors
            gamma.dh=gamma.dh, # from hospitalization to death
            theta1=theta1, # proportion infectious in the hospital
            delta1=delta1, # CFR for unhospitalized
            delta2=delta2, # CFR for hospitalize
            beta.I=beta.I * (1-z), # transmission rate in the community
            beta.H=beta.H * (1-z.H), # transmission rate in the hospital
            beta.F=beta.F * (1-z.F) # transmission rate at funerals, 
);

## simulation starts on Mar 1, no intervention until May 4 (~9 weeks later)
## last case identified on July 12 (~19w weeks later)
times1=seq(1,9*7,by=1); # first 9 weeks with out intervention
times2=seq(9*7, 20*7,by=1); # afterward with incidence

# Stage 1: no intervention
simNoCtrl=ode(y=state,times=times1,func=EBOLA,parms=parmsNoCtrl);

# Stage 2: yes intervention
state2=c(S=tail(simNoCtrl[,'S'],1),E=tail(simNoCtrl[,'E'],1),
         I=tail(simNoCtrl[,'I'],1),H=tail(simNoCtrl[,'H'],1),
         FF=tail(simNoCtrl[,'FF'],1),R=tail(simNoCtrl[,'R'],1),
         cumInci=tail(simNoCtrl[,'cumInci'],1));
simCtrl=ode(y=state2,times=times2,func=EBOLA,parms=parmsCtrl);

# Combine the two
sim2=rbind(simNoCtrl,simCtrl[-1,])

inci2=sim2[seq(7,nrow(sim2),by=7),'cumInci']-c(0,sim2[seq(7,nrow(sim2)-7,by=7),'cumInci']) # get weekly incidence


#######################################
## COMPARE: NO-INTERVENTION, YES-INTERVENTION, OBSERVATIONS
#######################################
## mock data
dates=seq(as.Date('1995/3/4'),as.Date('1995/7/15'),by=7);
obs=c(4,2,5,1,7,6,18,25,58,40,52,28,16,19,4,5,0,0,1,0); # 291 cases
da=cbind(dates,obs)

## plot results：simulated (with intervention) v. observations
par(mfrow=c(1,1), mar=c(3,3,1,1),mgp=c(1.6,.5,0),cex=1.1)
## Note: In barplot, bars are not drawn at intervals 1:10, etc., but on something else. 
## The actual x-coodinates can be seen if you save your barplot object. 
df.bar=barplot(obs,col='white',xlab='Time from onset (week)',
               ylab='Weekly incidence',ylim=c(0,80))
axis(1,at=df.bar,lab=1:20)
lines(x=df.bar, y=inci2,lwd=2) # use 'x=df.bar' to match with weeks
legend('topleft',c('Observed','Simulated: with intervention'),
       lty=c(NA,1),pch=c(0,NA),col=c('black','black'),cex=.8,bty='n')


##########################
### Q3: compare intervention with no intervention

  # CODE ON YOUR OWN plot and compare all three 





##############################################################################
## SENSITIVITY TESTS & EFFECTIVENESS OF INTERVENTION
##############################################################################
## Reduction in transmission within the community:
times1=seq(1,9*7,by=1); # first wk weeks with out intervention
times2=seq(9*7, 20*7,by=1);

zs=seq(.2,.8,by=.2)
res=matrix(0,20,length(zs))
for (i in 1:length(zs)){
  z=zs[i]; # intervention in the community
  z.H=0.75; # intervention in the hospital
  z.F=0.875; # intervention in safe burial
  parmsCtrl=c(alpha=alpha, # incubation period: 7 days
              gamma.h=gamma.h, # from onset to hopspitalization: 5 days
              gamma.d=gamma.d, # from onset to death: 9.6 days
              gamma.i=gamma.i, # from onset to end of infectiousness for survivors: 10 days
              gamma.f=gamma.f, # from death to traditional burial 2 days
              gamma.ih=gamma.ih, # from hospitalization to end of infectiousness for survivors
              gamma.dh=gamma.dh, # from hospitalization to death
              theta1=theta1, # proportion infectious in the hospital
              delta1=delta1, # CFR for unhospitalized
              delta2=delta2, # CFR for hospitalize
              beta.I=beta.I * (1-z), # transmission rate in the community
              beta.H=beta.H * (1-z.H), # transmission rate in the hospital
              beta.F=beta.F * (1-z.F) # transmission rate at funerals, 
  );
  simNoCtrl=ode(y=state,times=times1,func=EBOLA,parms=parmsNoCtrl);
  state2=c(S=tail(simNoCtrl[,'S'],1),E=tail(simNoCtrl[,'E'],1),I=tail(simNoCtrl[,'I'],1),
           H=tail(simNoCtrl[,'H'],1),FF=tail(simNoCtrl[,'FF'],1),R=tail(simNoCtrl[,'R'],1),
           Inci=tail(simNoCtrl[,'cumInci'],1));
  simCtrl=ode(y=state2,times=times2,func=EBOLA,parms=parmsCtrl);
  sim=rbind(simNoCtrl,simCtrl[-1,])
  # save the result before exiting the loop
  res[,i]=sim[seq(7,nrow(sim),by=7),'cumInci']-c(0,sim[seq(7,nrow(sim)-7,by=7),'cumInci']) # get weekly incidence
}

ymax=max(res)*1.2
par(mar=c(3,3,1,1),mgp=c(2,.5,0),cex=1.1)
df.bar=barplot(obs,col='white',xlab='Date of onset (week)',
               ylab='Weekly incidence',ylim=c(0,ymax))
axis(1,at=df.bar,lab=1:20)
matlines(x=df.bar, y=res,lty=1, lwd=2,col=rainbow(ncol(res))) # use 'x=df.bar' to match with weeks
legend(0,ymax/1.2,c(paste('z=',zs,sep='')),
       lty=1, lwd=2,col=rainbow(ncol(res)),bty='n',cex=.8)



## NOW ADAPT THE CODE ABOVE TO COMPLETE THE EXERCISE BELOW:


##############################
# intervention at the hospital
##############################


##############################
# Intervention at the funerals
##############################


##############################
## Timing of intervention
##############################
