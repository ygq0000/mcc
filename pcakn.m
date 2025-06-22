function [n_normal]=pcakn(data_normal, k)
     p = data_normal';
     %k = 10;
     m = size(p,2);
     n = zeros(3,m);
    %neighbors = k_nearest_neighbors(p, p, k+1);%  %matlab版本小于等于7.5用的
     neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1)); %matlab版本大于等于7.5用的
     for i = 1:m
    x = p(:,neighbors(2:end, i));
    p_bar = 1/k * sum(x,2);
    
    P = (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %spd matrix P
   
    [V,D] = eig(P);
    [~, idx] = min(diag(D)); % choses the smallest eigenvalue
    n_normal(:,i) = V(:,idx);   % returns the corresponding eigenvector    
     end
     n_normal=n_normal';
end