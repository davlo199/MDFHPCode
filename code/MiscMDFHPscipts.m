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

%% Compute mean intensity
%Need to load the correct estimate file (output of "MDFHIntensityNewSum.m")
tvec=Events;
Intensity=MDFHazardRate(Nband,tvec,DiscEvents,renewal,finalP,M0);

%%
mean((Intensity(1,ismember(tvec,DiscEvents(1,2:DiscEvents(1,1)+1))))) %Mean intensity of subprocess 1 at event times in subprocess 1
mean((Intensity(2,ismember(tvec,DiscEvents(2,2:DiscEvents(2,1)+1))))) %Mean intensity of subprocess 2 at event times in subprocess 2


%% Auxilliary functions
function lambda0=BackLamb(Nband,t,DiscEvents,renewal,Parameters)
lambda0=zeros(Nband,length(t));

        for k=1:Nband
        lambda0(k,:)=Parameters(k,1).*ones(1,length(t));
        

        end
end

function ExcitationMatrix=ExcitedIntensity(Parameters,DiscEvents,M0,t,Nband,renewal)
if renewal==1||renewal==2
    ADD=2;
else
    ADD=1;
end
ExcitationMatrix=zeros(Nband,Nband*length(t));
%Events can be some general time but only need matrix at event times, 
for I=1:Nband
    for J=1:Nband
        alpha=Parameters(I,J+ADD);
        gamma=Parameters(I,J+ADD+Nband);
        beta=Parameters(I,J+ADD+2*Nband);
        c=Parameters(I,J+ADD+3*Nband);
        D=DiscEvents(J,2:(DiscEvents(J,1)+1));
        MAG=DiscEvents(J+Nband,2:(DiscEvents(J,1)+1));
        if alpha==0
            continue
        else
        ExcitationMatrix(I,((J-1)*length(t)+1):J*length(t))=alpha.*FracExp(D,t,MAG,beta,gamma,c,M0);
        end
    
    end
end
1;
end

function Excited=FracExp(D,t,MAG,beta,gamma,c,M0)
%t is time vector to calculate matrix at
%D are the events of the band, MAG their respective magnitudes, M0 cutoff
%Events is time of all in catalogue
%magnitude
%beta and gamma parameters for entries
if length(D)~=length(MAG)
error('MAG corresponds to events in first argument D')
end
expMag=exp(gamma.*(MAG-M0));
MAT=zeros(length(D),length(t));
MAT2=MAT;
for j=1:length(t)
sdx=sum(D<t(j)); %j+1 column is jth event as forst is total # of events
if sdx==0
continue
else
MAT(1:sdx,j)=t(j)-D(1:sdx); 
end
end
IDX2=MAT~=0;
TempVEC=(c^beta).*ml(-(c.*MAT(IDX2)).^beta,beta,beta,1).*(MAT(IDX2).^(beta-1)).';
MAT2(IDX2)=TempVEC.';
Excited = expMag*MAT2;

end
function [Hazard,SumIntensity]=MDFHazardRate(Nband,t,DiscEvents,renewal,Parameters,M0)
%Calculates MDFHazrdrate at vector time t
AA=ExcitedIntensity(Parameters,DiscEvents,M0,t,Nband,renewal);
SumIntensity=zeros(Nband,length(t));
for JJ=1:Nband
SumIntensity=SumIntensity+AA(:,(JJ-1)*length(t)+1:JJ*length(t));
end
Hazard=BackLamb(Nband,t,DiscEvents,renewal,Parameters)+SumIntensity;

end