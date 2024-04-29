function [result] = l2p_norm(X,p)
%l_{2,p} norm 0<p<2
if p>0
    row_sum=(sum(X.^2,2)).^(p/2);
    result=(sum(row_sum));
end
if p==0
    row_sum=sum(X.^2,2);
    result=sum(row_sum~=0);
end
% result2= sum(sqrt(sum(X.^2,2)));% ||X||_21
end
