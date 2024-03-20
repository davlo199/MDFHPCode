%% Compute information gains
%Using the output from "PredCapETASMarked.m" and "PredCapMDFHP.m" below
%computes IGPT
Tmp=size(P);

X=zeros(NMagInt,NInt);
for I=1:NInt
T0=PredInt(I);
TF=PredInt(I+1);
if any(Events>T0&(Events<=TF))
    X(:,I)=logical(histcounts(MAG(Events>T0&(Events<=TF)),MAGInt));
end
end

Xlogic=logical(X);
%
bMag=1/mean(MAG-M0);
delI=diff(PredInt);
empiricalrate=length(Events)/Events(end);
PPois=zeros(size(P));
for i=1:length(PPois)
    for j=1:height(PPois)
        probinteg=(exp(-bMag*(MAGInt(j)-M0))-exp(-bMag*(MAGInt(j+1)-M0)))./(1-exp(-bMag*(MAGInt(end)-M0)));
        markedrate=empiricalrate*probinteg; %Rate of the refernece process for the kth magnitude band 
        PPois(j,i)=(1-exp(-empiricalrate*delI(i)))*probinteg; 
        %Probability of observing an event in the ith time interval in the
        %jth magnitude band for an iid exponential marked Poisson
    end
end

IG=zeros(2,NMagInt)';
for i1=1:NMagInt

ptmp=P(i1,:);
ppoistmp=PPois(i1,:);
XlogTmp=Xlogic(i1,:);

IG(i1,1)=sum(log(ptmp(XlogTmp)./ppoistmp(XlogTmp)));
IG(i1,2)=sum(log((1-ptmp(~XlogTmp))./(1-ppoistmp(~XlogTmp))));
end