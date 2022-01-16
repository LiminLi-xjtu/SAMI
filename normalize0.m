%X: p by n
%normalized each feature, such that each feature has zero mean and norm 1
function X = normalize0(X)

[p,n] = size(X);

l = sqrt(sum(X.*X,2));
id = find(l);
d = 1./l;
X(id,:) = X(id,:).*repmat(d(id),1,n);

end