function [vProductivity,mTransition] = tauchen_ram(na,delta,sigma,m)
% Purpose:    Finds a N-Markov chain that approximates an AR(1) process
%
%                y(t+1) = delta*y(t)+(1-delta^2)^(1/2)*eps(t+1)
%       where  eps(t+1) ~ N(0,sigma) 
% Format:     [a,mTransition] = tauchen_ram(N,delta,sigma,m)
%
% Input:      N       scalar, number of nodes for Y (or vProductivity), 
%             delta   scalar, = corr(y_t,y_{t+1})
%             sigma   scalar, std. dev. of epsilons
%             m       max +- std. devs.
%
% Output:     vProductivity    N*1 vector, nodes for productivity
%             mTransition      N*N matrix, transition probabilities
%
%    Rodrigo A. Morales Mendoza, nov 2019.
%
%    This procedure implements George Tauchen's algorithm,
%           following Joao Gomes' 937 Wharton lecture notes.

% Transition matrix

if na > 1
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
    vProductivity = vProductivity';
    mTransition   = P;
else
    vProductivity = 0;
    mTransition = 1;
end