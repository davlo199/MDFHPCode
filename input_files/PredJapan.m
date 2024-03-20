M0=4.75;
M=13;
Tprime=10;

CSV=1;

dataname='Japan';

MANx0=[]; %For manually using starting values

PAR=2; %PAR=0 for running on a nonparallel computer.
A=2600;
AA=2600;
B=4100;

NRep=2000;
renewal=0;
MagVec=[10,5.5,M0];
MAGInt=[M0,5.5,6.5,10];
Nband=2; %should be length(MagVec)-1

%MDFHP parameters
finalP =[    0.0287    0.0067    0.1120    2.1419    0.1187    0.7590    0.4245    5.4521    0.0333
    0.0974    0.0047    0.4722    2.7624    0.9593    0.8684    0.5314    2.4973    0.1519];

magparam =[    2.6360    2.6660];

%ETAS parameters
mu0=0.120;
Aprod=1.246;
del=1.597;
p=1.089;
cE=0.029;

