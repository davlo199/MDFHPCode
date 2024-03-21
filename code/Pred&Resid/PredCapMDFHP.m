function PredCapMDFHP(input_file)
run(input_file)
start=clock;
%%
%Reading in datafile
Data=strcat(dataname,'.csv');
Points=readtable(Data);
Time=Points(:,1);
Mag=Points(:,5);
MAG=table2array(Mag);
MAG=MAG(1:end);
TimeCell=table2array(Time);
TimeNum=cell2mat(TimeCell);
datestr=TimeNum(:,1:10);
timestr=TimeNum(:,12:(end-1));
for i=1:height(Time)
datestr(i,1:11)=pad(datestr(i,:),11);
end
dtstr=[datestr,timestr];
numdatetime=datenum(dtstr); %Absolute time of events in serial date which is number of days
Events=numdatetime; %Removes event at time 0, relative number of days since first event
idx=MAG>=M0;
Events=Events(idx).';
MAG=MAG(idx).';
Events=Events(A:B);
MAG=MAG(A:B);
Events=Events-Events(1);


DiscEvents=zeros(2*Nband,length(Events)+1);
for i=1:(Nband)
idx2=MAG>=MagVec(i+1)&MAG<MagVec(i);
DiscEvents(i,2:sum(idx2)+1)=Events(idx2);
DiscEvents(i+Nband,2:sum(idx2)+1)=MAG(idx2);
DiscEvents(i,1)=sum(idx2);
DiscEvents(i+Nband,1)=sum(idx2);
end

PredInt=[1e-6:2:Events(end), Events(end)]; %End points of time bins, each is 2 days

NInt=length(PredInt)-1; %Number of intervals for prediction

NMagInt=length(MAGInt)-1;


 %Number of simulations for each interal

P=zeros(NMagInt,NInt);
    for K=1:NInt
    P(:,K)=PredictionLoop(finalP, M0,Events, MAG, NRep, K, PredInt,Nband,MagVec,magparam,MAGInt).';
    if floor(K/100)==K/100
        disp(sprintf("Finished %d of %d",K,NInt))
    end
    end
runtime=clock-start;
basefilename=strcat(dataname,'MDFHPmarkPred','A',sprintf('%d',A),'B',sprintf('%d',B));

save(basefilename)

end

function OUT=PredictionLoop(finalP, M0,Events, MAG, NRep, kInt, PredInt,Nband,MagVec,magparam,MAGInt)


rng(kInt) %For reproducibility
T0=PredInt(kInt);
TF=PredInt(kInt+1);
    %History fot his interval
    TmpEvent=Events(Events<=T0);
    TmpMAG=MAG(Events<=T0);
     Prop=zeros(length(MAGInt)-1,NRep);

        for L=1:NRep
        magout=SimPrediction(finalP,M0,TmpEvent,TmpMAG,T0,TF,Nband,MagVec,magparam);
        Prop(:,L)=histcounts(magout,MAGInt);
        end

    Prop(Prop~=0)=1;
    OUT=sum(Prop,2)/NRep;


end


function OUT = SimPrediction(finalP,M0,Events,MAG,T0,TF,Nband,MagVec,magparam)
%Events and MAG constitute history, (T0,TF] is time interval we want to
%predict in
%Rest are parameters in the simulation
%Only need to determine whether there is one point in (T0,TF] since
%considering unmarked case.

renewal=0; %Background rate is Poisson. DO NOT USE any other settings
t = T0;
    
SimPoints = Events;
epsilon=1e-10;

DiscEvents=zeros(2*Nband,1);
OUT=[];
%
for i=1:(Nband)
idx2=MAG>=MagVec(i+1)&MAG<MagVec(i);
DiscEvents(i,2:sum(idx2)+1)=SimPoints(idx2);
DiscEvents(i+Nband,2:sum(idx2)+1)=MAG(idx2);
DiscEvents(i,1)=sum(idx2);
DiscEvents(i+Nband,1)=sum(idx2);
end

%
while t<=TF

    intTmp=MDFHazardRate(Nband,t+epsilon,DiscEvents,renewal,finalP,M0);
    M=sum(intTmp);
    E = exprnd(1 / M, 1, 1);
    t = t + E;
        if t>TF
            break
        end
    U = rand;
    intTmp2=MDFHazardRate(Nband,t,DiscEvents,renewal,finalP,M0);
    if(U<=sum(intTmp2)/M) %This is a point is accepted

       for ii1=1:Nband %This finds the lowest band the 
            if (U>sum(intTmp2(1:ii1-1))/M)&&(U<=sum(intTmp2(1:ii1))/M)
                break %Point is of type ii1
            end
       end


    DiscEvents(ii1,1)=DiscEvents(ii1,1)+1; %Adding an extra event to the counter
    DiscEvents(ii1+Nband,1)=DiscEvents(ii1+Nband,1)+1;


    DiscEvents(ii1,1+DiscEvents(ii1,1))=t;

    MAG=truncmag(magparam,MagVec,ii1);

    DiscEvents(ii1+Nband,1+DiscEvents(ii1,1))=MAG;
    OUT=[OUT,MAG]; %Only storing simulated magnitudes for proportions
    end

end
    
end

function OUT=truncmag(magparam,MagVec,ii1)

Bmag=magparam(ii1);
MMin=MagVec(ii1+1);
MMax=MagVec(ii1);

OUT=[];

while isempty(OUT)

    magtmp=exprnd(1/Bmag,1,1)+MMin;
    if magtmp<=MMax
        OUT=magtmp;
%         disp("Done")
    end
end
end



function lambda0=BackLamb(Nband,t,DiscEvents,renewal,Parameters)
lambda0=zeros(Nband,length(t));
if renewal==1
        for k=1:Nband
        x=[Parameters(k,1),Parameters(k,2)];
        N=DiscEvents(k,1);
        lambda0(k,:)=GamIntensity(x,DiscEvents(k,2:N+1),t);
        end
    elseif renewal==2
        for k=1:Nband
        x=[Parameters(k,1),Parameters(k,2)];
        N=DiscEvents(k,1);
        lambda0(k,:)=WBLIntensity(x,DiscEvents(k,2:N+1),t);
        end
    else
        for k=1:Nband
        lambda0(k,:)=Parameters(k,1).*ones(1,length(t));
        end

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
%TempVEC=(c^beta).*ml(-(c.*MAT(IDX2)).^beta,beta,beta,1).*(MAT(IDX2).^(beta-1)).';
ZMat2=MAT(IDX2); 
TempVEC=(c^beta).*(MLapp(ZMat2.',c,beta,beta).*(ZMat2.^(beta-1)).').';

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
