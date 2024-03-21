M0=4.75;
M=13;
Tprime=10;

CSV=1;
maxFval=2000;
dataname='Japan';
Nrand=200;
batchsize=1;
MANx0=[]; %For manually using starting values

PAR=2; %PAR=0 for running on a nonparallel computer.
A=2600;
B=4100;

NRep=1000;

renewal=0;

MagVec=[10,5.5,M0];
Nband=2; %should be length(MagVec)-1

parallel.batchSize=1;
parallel.inner=2;
parallel.outer=0;
parallel.slurmParameters={'--mem-per-cpu', '1000m', '--time', '30:00:00','--cpus-per-task','12'};
