% function [id,obj] = RAFG(X,alpha,beta,nClass,p,k)
function [W,id,obj] = RAFG(X,alpha,beta,nClass,p,k)
% Input
% X: dim*num data matrix
% alpha: coefficient of L2p
% beta: coefficient of adaptive graph
% nClass: number of clusters

maxiter = 20;
if nargin < 6
    k = 5;
end
    
if nargin < 5
    p = 1; % L_2p
end

num = size(X,2);
dim = size(X,1);

X0 = X';
mX0 = mean(X0);
X1 = X0 - ones(num,1)*mX0;
scal = 1./sqrt(sum(X1.*X1)+eps);
scalMat = sparse(diag(scal));
X = X1*scalMat;
X = X';

distX = L2_distance_1(X,X);
[distX1, idx] = sort(distX,2);
A = zeros(num);
rr = zeros(num,1);
for i = 1:num
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
infind=find(A>10);
A(infind)=1;
r = mean(rr);

A0 = (A+A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;

for iter = 1:maxiter
%     disp(iter);
    [W,F]=updateWF(X,L0,alpha,beta,num,dim,nClass,p);
    
    A = updateS(F,r,k,num);
    D = diag(sum(A));
    L0 = D-A;
   
    obj(iter) = l2p_norm(X'*W - F,1)+alpha*l2p_norm(W,p)+beta*(2*trace(F'*L0*F)+r*norm(A,'fro')^2);
end

sqW = (W.^2);
sumW = sum(sqW,2);
[~,id] = sort(sumW,'descend');
end
