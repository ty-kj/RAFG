function [W,F] = updateWF(X,L,alpha,beta,num,dim,nClass,p)
%L: laplacian matrix
%X: data matrix(dim*num)
%gamma: coefficient of L2p
%dim: dimension of X
%beta: coefficient of adaptive graph

INTER = 30;%100
G = eye(dim);
Q = eye(num);
epsilon=10^(-6);
for i = 1:INTER
    invQ = inv(Q+2*beta*L);

    M = X*(2*beta*Q*invQ*L)*X';
    tempXLXG = (M+alpha*G);
%     [vec,val] = eig(tempXLXG);
%     [~,di] = sort(diag(val));
%     W = vec(:,di(1:nClass));
    [W, ~, evs]=eig1(tempXLXG, nClass+1, 0);
    F = invQ*Q*X'*W;
        
    E = X'*W-F;
    tempQ = 0.5*(sum(E.^2,2)+epsilon).^(-0.5);
    Q = diag(tempQ);
    
    tempG = p/2 * (sum(W.^2,2)+epsilon).^(p/2-1);
    G = diag(tempG);

    w1(i) = trace(W'*M*W); % log Tr(WMW)
    w2(i) = alpha*l2p_norm(W,p);% gama*||W||_2p
    WResult(i) = w1(i)+w2(i);

%     if i > 1 && abs(WResult(i-1)-WResult(i)) < 0.000001 %1e-06
%         break;
%     end
    
    obj1(i) = l2p_norm(X'*W - F,1)+alpha*l2p_norm(W,p)+beta*2*trace(F'*L*F);
   
end
% plot(obj1(1:20));
% xlabel('Number of Iterations')
% ylabel('Objective Function')
end
