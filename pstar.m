function ps = pstar(pi)
%  solves for unconditional distribution 
%  for an NxN  Markov chain 

N   = length(pi);
eN  = zeros(1,N-1);
eN  =[eN 1];
J   = -ones(1,N-1);
IJ  = [eye(N-1) J'];
tpi = pi(:,1:(N-1));
%ps  = (-eN*tpi)*inv(IJ*tpi - eye(N-1));
ps  = (-eN*tpi)/(IJ*tpi - eye(N-1));
ps  = [ps (1 - ps*ones(N-1,1))];