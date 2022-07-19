cfg = coder.config("mex");
cfg.OptimizeReductions=1;
cfg.EnableAutoParallelization=1;
%cfg.EnableMexProfiling=1;
cfg.EnableImplicitExpansion=1;
%cfg.ExtrinsicCalls=true;
%cfg.EnableAutoParallelizationReporting=1;
Q=diag([1000,10,0,10,10,1]);
R=diag([0.01,0.01]);
Qn=diag([1000,10,0,10,10,1])*10;
Twip_main(Q,R,Qn,1)
% or just Twip_main()
pause
Disp("press enter to clear the Command Window and compile")
pause
codegen Twip_main -config cfg -args {Q,R,Qn,1}
Twip_main_mex(Q,R,Qn,1)