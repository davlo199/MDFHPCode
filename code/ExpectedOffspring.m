%Example for computing expected number of offspring for an event in
%subprocess i from a triggering event in subprocess j (Equation (9)).
%This example is the expected number of offspring in subprocess 1 from an
%event in subprocess 2 (i.e. i=2, j=1).

gam=3.583; %gamma_{12}
alpha=0.042; %alpha_{12}
B=7.839; %B_{22} for triggering event being in subprocess 1 use B_{11}
M0=4; 
M1=4.35;

ExpNumb=alpha*(exp((M1-M0)*(gam-B))-1)./((gam-B)*(1-exp(-B*(M1-M0)))); %Equation 9