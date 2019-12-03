function [vProductivity,mTransition] = rouwenhorstMatrix_ram(N,p,q,phi)
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
%   following http://www.econ.nyu.edu/user/violante/NYUTeaching/Macrotheory/Spring14/LectureNotes/lecture6_14.pdf

if N > 1
    
    [vProd,mTemp] = rouwenhorstMatrix_ram(N,p,q,phi);
    
    if N == 2
        [vProductivity,mTransition] = [p 1-p; q 1-q];
    else
        [vProd,mTemp] = rouwenhorstMatrix_ram(N-1,p,q,phi);
        
    end
    vProductivity = zeros(na,1);
    vProductivity(na)  = m*sigma;
    vProductivity(1)  = -vProductivity(na);
    ystep = (vProductivity(na)-vProductivity(1))/(na-1);
    for ii=2:(na-1)
        vProductivity(ii) = vProductivity(ii-1)+ystep;
    end 
    P   = eye(na);
    a_1 = vProductivity(1);
    a_n = vProductivity(na);
    for j = 1:na
        aj = vProductivity(j);
        upperBoundA = (a_1 - delta*aj + ystep/2)/(sigma*sqrt(1-delta^2));
        P(1,j) = normcdf(upperBoundA);
        lowerboundA = (a_n - delta*aj - ystep/2)/(sigma*sqrt(1-delta^2));
        P(na,j) = 1-normcdf(lowerboundA);
    end
    for i = 2:(na-1)
        for j = 1:(na)
            ai = vProductivity(i);
            aj = vProductivity(j);
            upperBoundA = (ai - delta*aj + ystep/2)/(sigma*sqrt(1-delta^2));
            lowerboundA = (ai - delta*aj - ystep/2)/(sigma*sqrt(1-delta^2));
            P(i,j) = normcdf(upperBoundA)-normcdf(lowerboundA);
        end
    end
    %P'
    mTransition   = P;
else
    vProductivity = 0;
    mTransition = 1;
end