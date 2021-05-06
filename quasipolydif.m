function D = quasipolydif(x,malpha,alphafunc)

%  The function DM =  poldif(x, maplha, B) computes the
%  differentiation matrices D1, D2, ..., DM on arbitrary nodes.

%  Input (non-constant weight):
%
%  x:        Vector of N distinct nodes.
%  malpha:   Vector of weight values alpha(x), evaluated at x = x(k).
%  B:        Matrix of size M x N,  where M is the highest 
%            derivative required.  It should contain the quantities 
%            B(ell,j) = beta(ell,j) = (ell-th derivative
%            of alpha(x))/alpha(x),   evaluated at x = x(j).

%  Output:
%     D:     DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.

%  J.A.C. Weideman, S.C. Reddy 1998

n = length(x);                      
       


alpha = malpha(:);                % Make sure alpha is a column vector
M = length(B(:,1));           % First dimension of B is the number 
                              % of derivative matrices to be computed
 
        I = eye(n);                  % Identity matrix.
        L = logical(I);              % Logical identity matrix.

       XX = x(:,ones(1,n));
       DX = XX-XX';                  % DX contains entries x(k)-x(j).

    DX(L) = ones(n,1);               % Put 1's one the main diagonal.

        c = alpha.*prod(DX,2);       % Quantities c(j).

        C = c(:,ones(1,n)); 
        C = C./C';                   % Matrix with entries c(k)/c(j).
   
        Z = 1./DX;                   % Z contains entries 1/(x(k)-x(j))
     Z(L) = zeros(n,1);              % with zeros on the diagonal.

        X = Z';                      % X is same as Z', but with 
     X(L) = [];                      % diagonal entries removed.
        X = reshape(X,n-1,n);

        Y = ones(n-1,n);             % Initialize Y and D matrices.
        D = eye(n);                  % Y is matrix of cumulative sums,
                                     % D differentiation matrices.
for ell = 1:M
        Y   = cumsum([B(ell,:); ell*Y(1:n-1,:).*X]); % Diagonals
        D   = ell*Z.*(C.*repmat(diag(D),1,n) - D);   % Off-diagonals
     D(L)   = Y(n,:);                                % Correct the diagonal
DM(:,:,ell) = D;                                     % Store the current D
end
