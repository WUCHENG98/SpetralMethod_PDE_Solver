function polyprime = basisdiff(poly)
% compute the coefficient matrix of {P_n'(x)}

% Input:
%    poly:  the coefficient matrix of the basis polynomial set {P_n(x)}

% Output:
%    polyprime: coefficient matrix of {P_n'(x)}

format long e

% find number of polynomials
[N, ~]= size(poly);

% find diffmatrix
A = zeros([1,N-1])';
B = diag(1:1:N-1);
C = 0;
D = zeros([1,N-1]);
diffmatrix = cell2mat({A,B;C,D});

% find polyprime
polyprime = diffmatrix*poly;
end