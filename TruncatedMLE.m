run(input_file)
%
%
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


Bmag=zeros(1,2);
negLlhoodout=Bmag;
SE=zeros(1,2);
for i=1:2
MAG2=DiscEvents(i+Nband,2:DiscEvents(i,1)+1);
M0=MagVec(i+1);
MMax=MagVec(i);
Nrand=1;
MLLik=zeros(4,Nrand).';

x0=2*rand()-1;
[MLE,Llhood,~,~,~,hessian]=fminunc(@(x)BMAG(M0,MAG2,x,MMax),x0); 
SE(i)=sqrt(inv(hessian));
Bmag(i)=exp(MLE);
negLlhoodout(i)=Llhood;
end


function OUT=BMAG(M0,MAG,x,MMax)
B=exp(x(1));
N=length(MAG);
OUT=-(N*log(B)-B*sum(MAG-M0)-N*log(1-exp(-B*(MMax-M0))));

end

