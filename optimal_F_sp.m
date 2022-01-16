function P=optimal_F_sp(P0,Y,X,F12, G)
opts.record = 0; %
opts.mxitr  = 1000;
opts.xtol = 1e-10;
opts.gtol = 1e-10;
opts.ftol = 1e-10;

% P0=orth(P0);
% tic; 
[P, out]= OptStiefelGBB(P0, @fun_F, opts, Y,X,F12,G); 
% tsolve = toc;
% out.fval = -2*out.fval; % convert the function value to the sum of eigenvalues
% fprintf('\nOptM: obj: %7.6e, itr: %d, nfe: %d, cpu: %f, norm(XT*X-I): %3.2e \n', ...
%             out.fval, out.itr, out.nfe, tsolve, norm(X'*X - eye(k), 'fro') );
end

% function [F, G] = fun(X,  A)
%   G = -(A*X);
%   F = 0.5*sum(dot(G,X,1));
% end

function [F, G] = fun_F(P ,Y,X,F12,G)
  G = 2*Y'*Y*P - 2*Y'*X + P*G;
  F = (norm(Y*P-X,'fro'))^2+ F12;
end

