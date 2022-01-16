% function [F, G] = fun(X,  A)
%   G = -(A*X);
%   F = 0.5*sum(dot(G,X,1));
% end

function [F, G] = fun(X ,W,Y)
  G = 2*X*W*W'-2*Y*W';
  F = (norm(X*W-Y,'fro'))^2;
end

