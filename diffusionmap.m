
function Y = diffusionmap(X,d)

sigma=200;
[T,~] = transition_matrix(X','classic',sigma);
[phi, ~] = eig_decompose_normalized(T,10);
Y=phi(:,2:d+1);


end


