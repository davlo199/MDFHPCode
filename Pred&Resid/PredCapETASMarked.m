function PredCapETASMarked(input_file)
run(input_file)
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

%Setting up simulation
%
dT=2; %Intervals of two day duration
PredInt=[1e-6:dT:Events(end), Events(end)]; %End points of time bins, each is 2 days


NInt=length(PredInt)-1; %Number of intervals fr prediction

NMagInt=length(MAGInt)-1;

NRep=2000; %Number of simulations for each interal

bMag=1/mean(MAG-M0); %MLE of iid exponential

P=zeros(NMagInt,NInt);
    for K=1:NInt
    P(:,K)=PredictionLoop(mu0, Aprod, del, p, M0, cE, Events, MAG, NRep, K, PredInt,bMag,MAGInt);
    if floor(K/50)==K/50
        disp(K)
    end
    end

basefilename=strcat(dataname,'ETASMarkedPrediction','A',sprintf('%d',A),'B',sprintf('%d',B));

save(basefilename)



end



function OUT=PredictionLoop(mu0, A, del, p, M0, cE, Events, MAG, NRep, kInt, PredInt,bMag,MAGInt)

    %rng(kInt) %For reproducibility
    T0=PredInt(kInt);
    TF=PredInt(kInt+1);
    %History for this interval
    TmpEvent=Events(Events<=T0);
    TmpMAG=MAG(Events<=T0);
    Prop=zeros(length(MAGInt)-1,NRep);
        for L=1:NRep
        [~,magout]=SimPrediction(mu0, A, del, p, M0, cE,TmpEvent,TmpMAG,T0,TF,bMag,MAGInt);
        Prop(:,L)=histcounts(magout,MAGInt);
        end
        Prop(Prop~=0)=1;
        OUT=sum(Prop,2)/NRep;

end


function [SimTimeOut,magout] = SimPrediction(mu0, A, del, p, M0, cE,Events,MAG,T0,TF,bMag,MAGInt)
%Events and MAG constitute history, (T0,TF] is time interval we want to
%predict in
%Rest are parameters in the simulation

    t = T0;
    epsilon = 1e-10;
    SimPoints = Events;
   
    % Simulate Points
   SimTimeOut=[];
   magout=[];
    while t<=TF
        M = mu0 + A * exp(del * (MAG - M0))*((1+(t+epsilon-SimPoints)/cE).^(-p)).' ;
        E = exprnd(1 / M, 1, 1);
        t = t + E;
        if t>TF
            break
        end
        U = rand;

        if (U < (mu0 + A * exp(del * (MAG - M0)) * ((1+(t-SimPoints)/cE).^(-p)).' / M))

            SimPoints=[SimPoints,t];
            magtmp=TruncMag(bMag,MAGInt);
            MAG=[MAG,magtmp];
            
            SimTimeOut=[SimTimeOut,t];

            magout=[magout,magtmp];
            
        end

    end
    
end

function OUT=TruncMag(bMag,MAGInt)
OUT=[];
while isempty(OUT)
    mag=exprnd(1/bMag,1,1)+MAGInt(1);
    if mag>MAGInt(end)
        continue
    else
        OUT=mag;
    end
end


end