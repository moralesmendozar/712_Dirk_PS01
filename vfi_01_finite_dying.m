function [VT,IGT,ST,CT] = vfi_01_finite_dying(k,A,QQ,bbeta,phi,r,...
    u,dinit,tol,maxiter,Vinit,T,ybar)
% VFI, solves infinity Bewley model for precautionary savings.
%       V - Value function               V(k,y)
%       IG - Index of policy function
%       S - savings policy function      s(k,y)
%       C - consumption policy function  c(k,y)
% Inputs:
%      phi contains probas of dying next period.
% Rodrigo Morales
%   November 2019.

nk      = length(k);
na      = length(A);

% initialize
K       = repmat(k,1,nk);  
KP      = repmat(k',nk,1);  % k prime 


% Record the V_t, IG_t, S_t and C_t for each iteration of time...
%   time = 1:60
VT = zeros(nk,na,T);
IGT = zeros(nk,na,T);
ST = zeros(nk,na,T);
CT = zeros(nk,na,T);

% note: TAKE ADVANTAGE OF MONOTONICTY AND CONVEXITY...

% value of dying is very negative, so better consume all today:
%Vtp1 = -1000*ones(nk,na); % represents the V_{t+1}, while V is the V_t
Vtp1 = zeros(nk,na); % represents the V_{t+1}, while V is the V_t

for tt = T:-1:1
    bbetaj = bbeta*phi(tt);
    if (mod(tt,10)==0 || tt ==1)
            fprintf(' Time iteration = %d \n', tt); 
    end
    
    % F contains the utility for all combinations (k,k,s)
    F       = zeros(nk,nk,na);

    for j = 1:na

        constemp = A(j)*ybar(tt) + (1+r)*K - KP;
        matUtilityAux = u( constemp );
        matUtilityAux(constemp<=0.0001) = -1e200;
        F(:,:,j) =matUtilityAux;

    end
    
    %looking for V, which is V_t, given V_{t+1}, which is VT
    V = Vinit;
    Vold    = zeros(nk,na);     % initial V for the loop to find V
    IG      = zeros(nk,na);     % index of policy
    iter    = 0;
    d = dinit;
    
    while d > tol && iter < maxiter 
        iter = iter + 1;
        for nProductivity = 1:na
            
            expectedValueFunction = Vtp1*QQ';
            
            % We start from previous choice (monotonicity of policy function)
            gridCapitalNextPeriod = 1;
            for nCapital = 1:nk
                %[V(i,j),IG(i,j)] = max((1-bbeta)*F(i,:,j)' + bbeta*Vold*QQ(j,:)');
                valueHighSoFar = -100000.0;
                nCapitalChoice = 1;
                %capitalChoice  = vGridCapital(nCapitalChoice);

                for nCapitalNextPeriod = gridCapitalNextPeriod:nk

                    valueProvisional = (1-bbetaj)*F(nCapital,nCapitalNextPeriod,nProductivity)+bbetaj*expectedValueFunction(nCapitalNextPeriod,nProductivity);

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        nCapitalChoice = nCapitalNextPeriod;
                        %capitalChoice = vGridCapital(nCapitalNextPeriod);
                        gridCapitalNextPeriod = nCapitalNextPeriod;
                    else
                        break; % We break when we have achieved the max
                    end    

                end

                V(nCapital,nProductivity) = valueHighSoFar;
                IG(nCapital,nProductivity) = nCapitalChoice; % Policy function
                %S(nCapital,nProductivity) = capitalChoice
            end
        end
        d = abs(max(max(V-Vold)));
        Vold = V;
        if (mod(iter,20)==0 || iter ==1)
            fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, d); 
        end
    end
    % policy
    S       = k(IG); 
    % Consumption:
    C = repmat(A,nk,1) + (1+r).*repmat(k,1,na) -S;
    
    VT(:,:,tt) = V;
    IGT(:,:,tt) = IG;
    ST(:,:,tt) = S;
    CT(:,:,tt) = C;
    Vtp1 = V;
    
end

