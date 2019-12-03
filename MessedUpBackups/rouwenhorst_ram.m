function [vProductivity,mTransition] = rouwenhorst_ram(N,p,q,Psi)
% Purpose:    Finds a N-Markov chain that approximates an AR(1) process
%           when the process is close tor be persistent (near unit root)
%
%     with      E[y(t+1)] = (q-p)Psi / (2-p-q)
%               V(y_t) = Psi(1-4s(1-s) + 4s(1-s)/(N-1) )
%               corr(yt+1,yt) = p + q - 1
%               E(y_t+1|y_k] = (q-p)Psi + (p + q - 1) y_k
% 
% Format:     [vProductivity,mTransition] = rouwenhorst_ram(N,p,q,phi)
%
% Approximate it through a discrete-space Markov chain over a symmetric
%       and evenly-spaced state space {y1 y2 ... yN} with:  -y1 = yN = Psi
%       Construct the Markov chain ThetaN in a recursive way.
%
% Input:      N     scalar, number of nodes for Y,
%             p     scalar, 
%             q     scalar,
%             Phi   scalar, limit Y_N
%
% Output:     vProductivity    N*1 vector, nodes for productivity
%             mTransition      N*N matrix, transition probabilities
%
%    Rodrigo A. Morales Mendoza, nov 2019.
%
%   following http://www.econ.nyu.edu/user/violante/NYUTeaching/Macrotheory/Spring14/LectureNotes/lecture6_14.pdf

if N > 1
    % Get the transition matrix:
    mTransition = rouwenhorstMatrix_ram(N,p,q);
    % get the vector or productivity shocks
    vProductivity = linspace(-Psi,Psi, N);
    
    %vProductivity = zeros(N,1);
%     vProductivity(N)  = Psi;
%     vProductivity(1)  = -vProductivity(N);
%     ystep = (vProductivity(N)-vProductivity(1))/(N-1);
%     for ii=2:(N-1)
%         vProductivity(ii) = vProductivity(ii-1)+ystep;
%     end
else
    vProductivity = 0;
    mTransition = 1;
end