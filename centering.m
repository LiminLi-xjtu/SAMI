%X: p by n
%normalized each feature, such that each feature has zero mean and norm 1
function X = centering(X)

[p,n] = size(X);
if center == 1
   xbar = mean(X,2);
   X = X - repmat(xbar,1,n);    
end

end