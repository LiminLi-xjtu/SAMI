function [F] = procrust2( A,B )
% PROCRUST Orthogonal Procrustes problem
%    A2 = PROCRUST( A, B ) applies an orthogonal transformation to matrix B
%    by multiplication with F such that A-B*F has minimum Frobenius norm. The
%    results B*F is returned as A2.
%
%    [A2, F] = PROCRUST( A, B ) also returns the orthogonal matrix F that was
%    used for the transformation.
%
%    [A2, F, R] = PROCRUST( A, B ) also returns the Frobenius norm of A-B*F.
% Author : E. Larsen
% Date   : 12/22/03
% Email  : erik.larsen@ieee.org
% Reference: Golub and van Loan, p. 601.
% Error checking
[m1, n1] = size(A);
[m2, n2] = size(B);
% msg = nargchk( 2, 2, nargin );
% if ~isempty( msg )
%     error( msg )
% end
% Do the computation
C = B.'*A;
[U, S, V] = svd( C );
F = V(:,1:n2)*U';


% A2 = A*F;
% Optional output of norm
% if nargout > 2
%     r = norm( B-A2, 'fro' );
% end
end
