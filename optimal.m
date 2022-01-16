function P=optimal(P0,W,X)
opts.record = 0; %
opts.mxitr  = 1000;
opts.xtol = 1e-10;
opts.gtol = 1e-10;
opts.ftol = 1e-10;

% P0=orth(P0);
% tic; 
[P, out]= OptStiefelGBB(P0, @fun, opts, W,X); 
% tsolve = toc;
% out.fval = -2*out.fval; % convert the function value to the sum of eigenvalues
% fprintf('\nOptM: obj: %7.6e, itr: %d, nfe: %d, cpu: %f, norm(XT*X-I): %3.2e \n', ...
%             out.fval, out.itr, out.nfe, tsolve, norm(X'*X - eye(k), 'fro') );
end