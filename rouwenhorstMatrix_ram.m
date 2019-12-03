function [mTransition] = rouwenhorstMatrix_ram(N,p,q)
%       Construct the Markov chain ThetaN matrix in a recursive way.
%    Approximate it through a discrete-space Markov chain over a symmetric
%       and evenly-spaced state space {y1 y2 ... yN} with:  -y1 = yN = Psi
%       Construct the Markov chain ThetaN in a recursive way.
% THIS FUNCTION SIMPLY GETS THE THETAN matrix RECURSIVELY :)
% Input:      N     scalar, number of nodes for Y,
%             p     scalar, 
%             q     scalar,
%             Phi   scalar, limit Y_N
%
%     with      E[y(t+1)] = (q-p)Psi / (2-p-q)
%               V(y_t) = Psi(1-4s(1-s) + 4s(1-s)/(N-1) )
%               corr(yt+1,yt) = p + q - 1
%               E(y_t+1|y_k] = (q-p)Psi + (p + q - 1) y_k
%
% Output:     mTransition      N*N matrix, transition probabilities
%
%    Rodrigo A. Morales Mendoza, nov 2019.
%
%   following http://www.econ.nyu.edu/user/violante/NYUTeaching/Macrotheory/Spring14/LectureNotes/lecture6_14.pdf

if N > 1
    if N == 2
        mTransition = [p, 1-p; 1-q, q];
    else
        mTemp = rouwenhorstMatrix_ram(N-1,p,q);
        M1 = [mTemp zeros(N-1,1); zeros(1,N)];
        M2 = [zeros(N-1,1) mTemp; zeros(1,N)];
        M3 = [zeros(1,N); mTemp, zeros(N-1,1)];
        M4 = [zeros(1,N) ; zeros(N-1,1), mTemp];
        mTransition = p*M1 + (1-p)*M2 + (1-q)*M3 + q*M4;
        mTransition(2:end-1,:) = mTransition(2:end-1,:)./2;
    end
else
    mTransition = 1;
end