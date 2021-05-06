function P = polyev(pts, poly)
% evaluate the polynomial sets at the given points

% Input:
%    pts:   quadrature pts/nodes
%    poly:  the coefficient matrix of the basis polynomial set
% note: length of pts should be the same as the width/length of poly

% Output:
%    P:     the matrix of values of polynomial at given points (row: function //column: point)

format long e

% dertermine number of basis
[N, ~]= size(poly);

% generate the matrix of powers
n = (0:1:N-1)';
power = zeros(N);
for i = 1:1:N
    power(:,i) = pts(i).^n;
end    

P = (poly'*power)';
end 