hold on
if renewal==1||renewal==2
    ADD=2;
else
    ADD=1;
end
COMP=zeros(Nband,length(Events));

for II=1:Nband
    
    for JJ=1:length(Events)
        
        BackInt=Events(JJ).*finalP(II,1);

        disp(JJ)
       
        COMP(II,JJ)=BackInt+ExciteInt(Nband,DiscEvents,finalP,M0,Events(JJ),renewal,II); 
    
    end
    
    
end

%%

EVDX=ismember(Events,DiscEvents(1,2:DiscEvents(1,1)+1));

mdtau1=COMP(1,EVDX); %Transformed inter-event time residual process for subprocess 1

EVDX=ismember(Events,DiscEvents(2,2:DiscEvents(2,1)+1));

mdtau2=COMP(2,EVDX); %Transformed inter-event time residual process for subprocess 2



%% Auxilliary functions
function lambda0=BackLamb(Nband,t,DiscEvents,renewal,Parameters)
lambda0=zeros(Nband,length(t));

        for k=1:Nband
        lambda0(k,:)=Parameters(k,1).*ones(1,length(t));
        end

end

function ExcitationMatrix=ExcitedIntensity(Parameters,DiscEvents,M0,t,Nband,renewal)

    ADD=1;

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

function ExcitedIntegral=ExciteInt(Nband,DiscEvents,Parameters,M0,TF,renewal,I)
ExcitedIntegral=0;

    ADD=1;


    for J=1:Nband
        alpha=Parameters(I,J+ADD);
        gamma=Parameters(I,J+ADD+Nband);
        beta=Parameters(I,J+ADD+2*Nband);
        c=Parameters(I,J+ADD+3*Nband);
        D=DiscEvents(J,2:(DiscEvents(J,1)+1));
        D=D(D<=TF);
        MAG=DiscEvents(J+Nband,2:(DiscEvents(J,1)+1));
        if alpha==0
            continue
        else
        betaVEC=(ml(-(c.*(TF-D)).^beta,beta)-1).';
        expVEC=exp(gamma*(MAG(1:length(D))-M0));
        ExcitedIntegral=ExcitedIntegral-alpha.*expVEC*betaVEC;
        end
    end

end
