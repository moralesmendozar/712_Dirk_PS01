function [V,IG,S,C] = vfi_01_infty(k,A,QQ,bbeta,r,u,d,tol,maxiter,V)
% VFI, solves infinity Bewley model for precautionary savings.
%       V - Value function               V(k,y)
%       IG - Index of policy function
%       S - savings policy function      s(k,y)
%       C - consumption policy function  c(k,y)
% Rodrigo Morales
%   November 2019.

nk      = length(k);
na      = length(A);

% initialize
K       = repmat(k,1,nk);  
KP      = repmat(k',nk,1);  % k prime 
Vold    = zeros(nk,na);     % initial V
IG      = zeros(nk,na);     % index of policy

% F contains the utility for all combinations (k,k,s)
F       = zeros(nk,nk,na);

for j = 1:na
    F(:,:,j) = u( A(j) + (1+r)*K - KP );
end


% VFI
% iter    = 0;
% while d > tol && iter < maxiter 
%     iter = iter + 1;
%     for i = 1:nk
%         for j = 1:na 
%             [V(i,j),IG(i,j)] = max((1-bbeta)*F(i,:,j)' + bbeta*Vold*QQ(j,:)');
%         end
%     end
%     d = abs(max(max(V-Vold)));
%     Vold = V;
%     if (mod(iter,10)==0 || iter ==1)
%         fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, d); 
%     end
% end

% TAKE ADVANTAGE OF MONOTONICYT AND CONVEXITY...

iter    = 0;
while d > tol && iter < maxiter 
    iter = iter + 1;
    for j = 1:na
        
        %CHECK ALL THISSSSSS 
        
        
        
        
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        for nCapital = 1:nk
            %[V(i,j),IG(i,j)] = max((1-bbeta)*F(i,:,j)' + bbeta*Vold*QQ(j,:)');
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nk
                
                %consumption = mOutput(nCapital,nProductivity)-vGridCapital(nCapitalNextPeriod);
                prodctvt = exp(vProductivity(nProductivity));
                currentCapital = vGridCapital(nCapital);
                kprimeTomorow = vGridCapital(nCapitalNextPeriod);
                invstm = investment(prodctvt,currentCapital,kprimeTomorow);
                %for 01.(d)(3)
                if invstm > delta*currentCapital
                    dividend = profitMatrix(nCapital,nProductivity) - invstm - phi(prodctvt,currentCapital,kprimeTomorow);
                else
                    dividend = profitMatrix(nCapital,nProductivity) - invstm - phi2(prodctvt,currentCapital,kprimeTomorow);
                end
                valueProvisional = (1-bbeta)*dividend+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
            
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    
                  
            end
            
            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;
            
        end
    end
    d = abs(max(max(V-Vold)));
    Vold = V;
    if (mod(iter,10)==0 || iter ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, d); 
    end
end


% policy
S       = k(IG); 

% Consumption:
C = repmat(A',nk,1) + (1+r).*repmat(k,1,na) -S;
