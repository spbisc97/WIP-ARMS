cfg = coder.config("mex");
cfg.OptimizeReductions=1;
cfg.EnableAutoParallelization=1;
%cfg.EnableMexProfiling=1;
cfg.EnableImplicitExpansion=1;
%cfg.ExtrinsicCalls=true;
%cfg.EnableAutoParallelizationReporting=1;
Q=diag([10,1,10,1,10,1]);
R=diag([0.01,0.01]);
Qn=diag([10,1,10,1,10,1])*1000;

codegen Twip_main -config cfg -args {Q,R,Qn,1}
Twip_main_mex(Q,R,Qn,1)