function D = polydif(pts,wts,poly)
% compute the derivative matrix with repesct to the physical space

% Input:
%    pts:   quadrature pts/nodes
%    wts:   quadrature wts
%    poly:  the coefficient matrix of the basis polynomial set

% Output:
%    D:     the derivative matrix

format long e

[N, ~]= size(poly);

% evaluate the coefficient matrix of {P_n'(x)}
polyprime = basisdiff(poly);

% evluate poly and polyprime at the given pts
P = polyev(pts, poly);
Pprime = polyev(pts, polyprime);

% evaluate the derivative matrix
D = zeros(N);
for i=1:N
    for j=1:N
        s=sum(P(j,:).*Pprime(i,:)); %sum over P'_n(x_i)P_n(x_j)
        D(i,j)=sqrt(wts(i)*wts(j))*s;
    end
end
end


