function [mTransition] = rouwenhorst_ram(N,p,q,phi)
% Purpose:    Finds a N-Markov chain that approximates an AR(1) process
%           when the process is close tor be persistent (near unit root)
%
%     with      E[y(t+1)] = (q-p)Psi / (2-p-q)
%               V(yt) = Psi(1-4s(1-s) + 4s(1-s)/(N-1) )
%       where  eps(t+1) ~ N(0,sigma) 
% Format:     [a,mTransition] = tauchen_ram(N,delta,sigma,m)
%
% Approximate it through a discrete-space Markov chain over a symmetric
%       and evenly-spaced state space {y1 y2 ... yN} wit -y1 = yN = Psi
%       Construct the Markov chain ThetaN in a recursive way.
%
% Input:      N       scalar, number of nodes for Y,
%             delta   scalar, corr(y_t,y_{t+1})
%             sigma   scalar, std. dev. of epsilons
%             m       max +- std. devs.
%
% Output:     vProductivity    N*1 vector, nodes for productivity
%             mTransition      N*N matrix, transition probabilities
%
%    Rodrigo A. Morales Mendoza, nov 2019.
%
%   following 

if N > 1
    if N == 2
        mTransition = [p 1-p; q 1-q];
    else
        mTemp = rouwenhorstMatrix_ram(N-1,p,q,phi);
        M1 = [mTemp zeros(N-1,1); zeros(1,N-1), 0];
        M2 = [zeros(N-1,1) mTemp; 0  zeros(1,N-1)];
        M3 = [zeros(1,N-1) 0; mTemp, zeros(N-1,1)];
        M4 = [0 zeros(1,N-1) ; zeros(N-1,1), mTemp];
        mTransition = p*M1 * (1-p)*M2 + q*M3 + (1-q)*M4;
    end
else
    mTransition = 1;
end