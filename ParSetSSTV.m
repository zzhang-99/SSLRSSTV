function  [opts]=ParSetSSTV
opts.mu= 0.1;
opts.lambda =0.3; 
opts.beta1 =1; 
opts.tau=0.03;
opts.r1 =2;opts.r2 =1;
opts.blocksize = 20;opts.stepsize  = 10; 
opts.maxIter = 50;opts.Innerloop_B = 1;
opts.tol = 1e-6;