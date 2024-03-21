M0=4;
M=13;
Tprime=10;

CSV=1;

dataname='MAT7623';


PAR=2; %PAR=0 for running on a nonparallel computer.
A=1259;
B=5393;

NRep=2000;

renewal=0;

MagVec=[10,4.35,M0];
MAGInt=[M0,4.35,5.35,10];
Nband=2; %should be length(MagVec)-1

%MDFHP parameters
finalP = [ 0.0798    0.0408    0.0423    1.3330    3.5828    0.7176    0.6684   11.4691    0.4618
    0.0488    0.1155    0.8075    1.2068    0.3916    0.6867    0.6232    2.5828    0.0654];

magparam=[2.4687    7.8385];

%ETAS parameters
mu0=0.119;
Aprod=1.767;
del=1.135;
p=0.962;
cE=0.022;
