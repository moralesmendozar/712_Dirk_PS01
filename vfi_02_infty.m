function [V,IG,S,C] = vfi_02_infty(k,A,QQ,bbeta,r,u,d,tol,maxiter,V)
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
    constemp = A(j) + (1+r)*K - KP;
    matUtilityAux = u( constemp );
    matUtilityAux(constemp<=0.00001) = -1e200;
    F(:,:,j) =matUtilityAux;
%     F(:,:,j) = u( constemp );
%     % if negative consumption, fix very large negative utility:
%     for nii = 1:nk
%         for njj = 1:nk
%             if constemp(nii,njj) <=0.0001
%                 F(nii,njj,j) = -1e-200;
%             end
%         end
%     end
end


% VFI
iter    = 0;
while d > tol && iter < maxiter 
    iter = iter + 1;
    for j = 1:na
        for i = 1:nk
            [V(i,j),IG(i,j)] = max((1-bbeta)*F(i,:,j)' + bbeta*Vold*QQ(j,:)');
        end
    end
    d = max(max(abs(V-Vold)));
    %d = norm(V-Vold,1);
    Vold = V;
    if (mod(iter,10)==0 || iter ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, d); 
    end
end

% TAKE ADVANTAGE OF MONOTONICTY AND CONVEXITY...

% iter    = 0;
% while d > tol && iter < maxiter 
%     iter = iter + 1;
%     for nProductivity = 1:na
%
%         expectedValueFunction = V*QQ;
%         
%         % We start from previous choice (monotonicity of policy function)
%         gridCapitalNextPeriod = 1;
%         for nCapital = 1:nk
%             %[V(i,j),IG(i,j)] = max((1-bbeta)*F(i,:,j)' + bbeta*Vold*QQ(j,:)');
%             valueHighSoFar = -1000.0;
%             nCapitalChoice = 1;
%             %capitalChoice  = vGridCapital(nCapitalChoice);
%             
%             for nCapitalNextPeriod = gridCapitalNextPeriod:nk
%                 
%                 valueProvisional = (1-bbeta)*F(nCapital,nCapitalNextPeriod,nProductivity)+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
%             
%                 if (valueProvisional>valueHighSoFar)
%                     valueHighSoFar = valueProvisional;
%                     nCapitalChoice = nCapitalNextPeriod;
%                     %capitalChoice = vGridCapital(nCapitalNextPeriod);
%                     gridCapitalNextPeriod = nCapitalNextPeriod;
%                 else
%                     break; % We break when we have achieved the max
%                 end    
%                   
%             end
%             
%             V(nCapital,nProductivity) = valueHighSoFar;
%             IG(nCapital,nProductivity) = nCapitalChoice; % Policy function
%             %S(nCapital,nProductivity) = capitalChoice
%         end
%     end
%     d = abs(max(max(V-Vold)));
%     Vold = V;
%     if (mod(iter,10)==0 || iter ==1)
%         fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, d); 
%     end
% end


% policy
S       = k(IG); 

% Consumption:
C = repmat(A,nk,1) + (1+r).*repmat(k,1,na) -S;
