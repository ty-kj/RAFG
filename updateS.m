function A1 = updateS(F,r,k,num)

% method-1:
% distf = L2_distance_1(F',F');
% distf = sqrt(distf);
% [distX1, idx] = sort(distf,2);
% A = zeros(num);
% for i=1:num
%     idxa0 = idx(i,2:k+1);
%     dfi = distf(i,idxa0);
%     
%     ad = -(dfi)/(2*r);
%     A(i,idxa0) = EProjSimplex_new(ad);
% end
% 
% infind = find(A>10);
% A(infind) = 1;
% 
% A = (A+A')/2;
% A1= 0.5*A./(distf+eps);

% A1=0.5*A.*(distf+eps).^(-1);
% tempQ = 0.5* (sqrt(sum(E.^2,2)+eps)).^(-1);

% method-2:
distf = L2_distance_1(F',F');
distf = sqrt(distf);
[~, idx] = sort(distf,2);
A = zeros(num);
for i=1:num
    idxa0 = idx(i,2:k+1);
    dfi = distf(i,idxa0);
    
    ad = -(dfi)/(2*r);
    A(i,idxa0) = EProjSimplex_new(ad);
end

infind = find(A>10);
A(infind) = 1;
A1 = (A+A')/2;
end
