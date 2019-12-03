function [lambda] = findLambdaStationary_MATRIX(lambdainit,matTransttn,IG,tol, maxiter)

d = 100;
[nk, ns] = size(IG);
lambda = lambdainit';
lambda = lambda(:); %vectorize the matrix
lambdainit = lambda;
Gamma = sparse(nk*ns,nk*ns);
iter = 0;

for jj = 1:nk
    for ii = 1:ns
        pis = matTransttn(ii,:);
        Gamma((jj-1)*ns + ii, (IG(jj,ii)-1)*ns+1:IG(jj,ii)*ns) = pis';
    end
end

while d > tol && iter<maxiter
    iter = iter+1;
    
    lambda = (Gamma')*lambdainit;
    d = norm(lambdainit-lambda,2);
    lambdainit = lambda;
    
%     if(mod(iter,20)==0)
%         fprintf('lambdaStat Iteration = %d, Sup Diff = %2.8f\n', iter, d); 
%     end
end

lambda = reshape(lambda,ns,nk);
lambda = lambda';

end