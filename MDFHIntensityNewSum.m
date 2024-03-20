function MDFHIntensityNewSum(input_file)
run(input_file)
%% Formulate Earthquake data and check data
Nband=length(MagVec)-1;

%First column is a, then b, then matrix of alpha_ij, beta_ij, gamma_ij as in model same order

if CSV==1
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
else
    Data=strcat(dataname,'.txt');
    Points=readtable(Data);
    Time=Points(:,2);
    Mag=Points(:,6);
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
end
start=datetime("now");
%
DiscEvents=zeros(2*Nband,length(Events)+1);
for i=1:(Nband)
idx2=MAG>=MagVec(i+1)&MAG<MagVec(i);
DiscEvents(i,2:sum(idx2)+1)=Events(idx2);
DiscEvents(i+Nband,2:sum(idx2)+1)=MAG(idx2);
DiscEvents(i,1)=sum(idx2);
DiscEvents(i+Nband,1)=sum(idx2);
end



%
if renewal==1||renewal==2
    ADD=2;
else
    ADD=1;
end

%When actually doing minimisation need to check that there are more
%data points than number of parameters for each magnitude band.
Nparam=ones(1,Nband);

Nparam=(ADD+4.*Nband).*Nparam;
%Check more data than parameters

for j=1:(Nband)
jdx=DiscEvents(j,:)~=0;
if sum(jdx)<Nparam(j)
    error('Less data than parameters in band [%d,%d] change band or data set',MagVec(j),MagVec(j+1))
end
end
totalParam=sum(Nparam);



% Slurm Parallelisation

if PAR==2
functionhandle=@(Nrand)(LOOPx0minimise(Nrand,batchsize,Nband,DiscEvents,renewal,M0,Events,Events(end),MANx0,maxFval));
finalParametersMatrix=slurmParallel(functionhandle,1:Nrand,parallel.slurmParameters);
MLLik=finalParametersMatrix.';
else
    %MLLik=zeros(Nrand,1+totalParam);
    for K=1:Nrand
    MLLik(K,:)=LOOPx0minimise(Nrand,batchsize,Nband,DiscEvents,renewal,M0,Events,Events(end),MANx0,maxFval);
    end

end
%




jdx=MLLik(:,end)==min(MLLik(:,end));
if sum(jdx)>1
AA=MLLik(jdx,:);
TransEst=AA(1:totalParam,:);
else 
TransEst=MLLik(jdx,:);
TransEst=TransEst(1:totalParam);
end

IDX=[true(Nband,ADD),false(Nband),true(Nband),false(Nband),true(Nband)];
finalP(IDX)=exp(TransEst(IDX));
finalP(~IDX)=exp(TransEst(~IDX))./(1+exp(TransEst(~IDX)));
finalP=reshape(finalP,Nband,ADD+4*Nband);
BestMinlik=min(MLLik(:,end));

basefilename=strcat(dataname,'MDFHP','Nband',sprintf('%d',Nband),'A',sprintf('%d',A),'B',sprintf('%d',B));

endtime=datetime("now");

runtime=endtime-start; 
save(basefilename)

end

function OUT=LOOPx0minimise(nRand,batchsize,Nband,DiscEvents,renewal,M0,Events,TF,MANx0,maxFval)
    for run=1:batchsize
        nRun=nRand*batchsize+run;
        rng(nRun)
        if isempty(MANx0)
        if renewal~=1||renewal~=2
        Par=zeros(Nband,1+4*Nband); 
        Par(1:Nband,1)=-3+6*rand(Nband,1);
        Par(1:Nband,2:Nband+1)=-3+6*rand(Nband);
        Par(1:Nband,Nband+2:2*Nband+1)=-3+6*rand(Nband);
        Par(1:Nband,2*Nband+2:3*Nband+1)=-3+6*rand(Nband);
        Par(1:Nband,3*Nband+2:end)=-3+6*rand(Nband);

        else
        Par=zeros(Nband,2+4*Nband); 
        Par(1:Nband,2)=-3+6*rand(Nband,1);
        Par(1:Nband,1)=-3+6*rand(Nband,1);
        Par(1:Nband,3:Nband+2)=-3+6*rand(Nband);
        Par(1:Nband,Nband+3:2*Nband+2)=-3+6*rand(Nband);
        Par(1:Nband,2*Nband+3:3*Nband+2)=-3+6*rand(Nband);
        Par(1:Nband,3*Nband+3:end)=-3+6*rand(Nband);
        end
        else

        Par=MANx0+0.1*randn(size(MANx0));

        end
        Parameters0=reshape(Par,1,numel(Par));
        [MLE,Llhood,SE]=MDFHPMinimiser(Parameters0,Nband,DiscEvents,renewal,M0,Events,TF,maxFval);
        SE=SE.';
        OUT=[MLE,SE,Llhood];       
    end
end


%Minimisation
function [MLE,Llhood,SE]=MDFHPMinimiser(Parameters0,Nband,DiscEvents,renewal,M0,Events,TF,maxFval)

options=optimoptions('fminunc','MaxFunctionEvaluations',maxFval,'Display','iter','StepTolerance',1e-6,'OptimalityTolerance',1e-6);

%[MLE,Llhood,~,~,~,hessian]=fminunc(@(Parameters)MDFHnegLogLik(Nband,DiscEvents,renewal,Parameters,M0,Events,TF),Parameters0,options); 
[MLE,Llhood,~,~,~,~]=fminunc(@(Parameters)MDFHnegLogLik(Nband,DiscEvents,renewal,Parameters,M0,Events,TF),Parameters0,options); 
hessian=eye(length(MLE));
SE=sqrt(diag(inv(hessian)));
end
%% Likelihood function
function neglogLike=MDFHnegLogLik(Nband,DiscEvents,renewal,Parameters,M0,Events,TF)
if renewal~=1&&renewal~=2 %Background rate is Poisson
Parameters=reshape(Parameters,[Nband,1+4*Nband]);
IDX=[true(Nband,1),false(Nband),true(Nband),false(Nband),true(Nband)];
TransPar(IDX)=exp(Parameters(IDX));
TransPar(~IDX)=exp(Parameters(~IDX))./(1+exp(Parameters(~IDX)));
TransPar=reshape(TransPar,Nband,1+4*Nband);
else
Parameters=reshape(Parameters,[Nband,2+4*Nband]);
IDX=[true(Nband,2),false(Nband),true(Nband),false(Nband),true(Nband)];
TransPar(IDX)=exp(Parameters(IDX));
TransPar(~IDX)=exp(Parameters(~IDX))/(1+exp(Parameters(~IDX)));
TransPar=reshape(TransPar,Nband,2+4*Nband);
end
logSUM=0;
for KK=1:Nband
logSUM=logSUM+sum(log(MDFHazardRate(Nband,DiscEvents(KK,2:DiscEvents(KK,1)+1),DiscEvents,renewal,TransPar,M0,KK)));
end



    BackInt=Events(end).*(sum(TransPar(:,1)));

CompInt=BackInt+ExciteInt(Nband,DiscEvents,TransPar,M0,TF,renewal); 
neglogLike=CompInt-logSUM;

end


%% Various auxillary functions for likelihood
function ExcitedIntegral=ExciteInt(Nband,DiscEvents,Parameters,M0,TF,renewal)
ExcitedIntegral=0;
if renewal==1||renewal==2
    ADD=2;
else
    ADD=1;
end
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
        betaVec=zeros(size(D));

        cTime=(c.*(TF-D)).^beta;

        Intdx1=cTime>=10; %Events to use Poin Asym App.

        betaVec(Intdx1)=AppML(13,beta,-cTime(Intdx1),1);

        if beta>=0.5

        Intdx2=cTime<=4;

        betaVec(Intdx2)=PSML(120,beta,-cTime(Intdx2),1);

        Intdx3=~(Intdx1|Intdx2);

        betaVec(Intdx3)=ml(-cTime(Intdx3),beta);    

        elseif beta<0.5 && beta>=0.25

        Intdx2=cTime<=2;

        betaVec(Intdx2)=PSML(200,beta,-cTime(Intdx2),1);

        Intdx3=~(Intdx1|Intdx2);

        betaVec(Intdx3)=ml(-cTime(Intdx3),beta);

        else

        Intdx2=cTime<=1;
    
        betaVec(Intdx2)=PSML(100,beta,-cTime(Intdx2),1);

        Intdx3=~(Intdx1|Intdx2);

        betaVec(Intdx3)=ml(-cTime(Intdx3),beta);
    
        end
        
        expVEC=exp(gamma*(MAG-M0));

        ExcitedIntegral=ExcitedIntegral-alpha.*expVEC*(betaVec-1).';
        end
    end
end
end
function Hazard=MDFHazardRate(Nband,t,DiscEvents,renewal,Parameters,M0,I)
%Calculates MDFHazrdrate at vector time t in I band

Hazard=BackLamb(t,DiscEvents,renewal,Parameters,I)+ExcitedIntensity(Parameters,DiscEvents,M0,t,Nband,renewal,I);

end



%% Exciting functions
function ExcitationMatrix=ExcitedIntensity(Parameters,DiscEvents,M0,t,Nband,renewal,I)

    ADD=1;

ExcitationMatrix=zeros(1,length(t));
%Events can be some general time but only need matrix at event times, 

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
        ExcitationMatrix=ExcitationMatrix+alpha.*FracExp(D,t,MAG,beta,gamma,c,M0);
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

ZMat=(c.*MAT(IDX2)).^beta;

ZMat2=MAT(IDX2); 

SumDx1=ZMat>=10;

TVec1=zeros(size(ZMat2));

TVec1(SumDx1)=(c^beta).*(AppML(13,beta,-ZMat(SumDx1).',beta).*(ZMat2(SumDx1).^(beta-1)).').';

if beta>=0.5
    
    SumDx2=ZMat<=4;

    TVec1(SumDx2)=(c^beta).*(PSML(120,beta,-ZMat(SumDx2).',beta).*(ZMat2(SumDx2).^(beta-1)).').';

    SumDx3=~(SumDx1|SumDx2);

    TVec1(SumDx3)=(c^beta).*ml(-ZMat(SumDx3),beta,beta).*(ZMat2(SumDx3).^(beta-1)).';

elseif beta<0.5 && beta>=0.25

    SumDx2=ZMat<=2;

    TVec1(SumDx2)=(c^beta).*(PSML(200,beta,-ZMat(SumDx2).',beta).*(ZMat2(SumDx2).^(beta-1)).').';

    SumDx3=~(SumDx1|SumDx2);

    TVec1(SumDx3)=(c^beta).*ml(-ZMat(SumDx3),beta,beta).*(ZMat2(SumDx3).^(beta-1)).';

else

    SumDx2=ZMat<=1;

    TVec1(SumDx2)=(c^beta).*(PSML(100,beta,-ZMat(SumDx2).',beta).*(ZMat2(SumDx2).^(beta-1)).').';

    SumDx3=~(SumDx1|SumDx2);

    TVec1(SumDx3)=(c^beta).*ml(-ZMat(SumDx3),beta,beta).*(ZMat2(SumDx3).^(beta-1)).';

end

MAT2(IDX2)=TVec1.';

Excited = expMag*MAT2;

end



function lambda0=BackLamb(t,DiscEvents,renewal,Parameters,k)


        
    lambda0=Parameters(k,1).*ones(1,length(t));
        

end




function OUT=AppML(M,beta1,Z,beta2)

% OUT=0;
% for II=1:M
% OUT=OUT-(Z.^(-II))./gamma(beta-II.*beta);
% end
%I think below is vectorised version of above

OUT=-(gamma(beta2-beta1.*(1:M)).^(-1))*(Z.^(-(1:M).'));

end

function OUT2=PSML(M,beta1,Z,beta2)
% OUT=0;
% for k=0:M
% OUT=OUT+(Z.^k)./gamma(beta2+k.*beta1);
% end
%I think below is vectorised version of above
OUT2=real(sum(exp((0:M).*log(Z).'-gammaln(beta2+beta1.*(0:M))),2)).';
OUT2(abs(Z)<=1e-15)=1./gamma(beta2);
end