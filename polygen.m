function poly = polygen(name,N)
% generate the coefficient matrix given the name of the .dat file

% Input:
%   name: the name of the .dat file containing the alphas and betas
%   N:    the number of basis polynomial required to be generated

% Output:
%   poly: The coefficient matrix of the basis poly set

% load alpha and beta (and mu0)
fileName= ['ab_',name,'.dat'];
data = importdata(fileName);
a = data(:,1);
b = data(:,2);
mu0 = data(1,3);
rtb = sqrt(b);
rtb(1)=[];

% initialize poly coeff matirx
poly = zeros(N);

% coeff of P_1
poly(1,1) = 1/sqrt(mu0);

% coeff of P_2
poly(1,2) = -a(1)/sqrt(b(2)*mu0);
poly(2,2) = 1/sqrt(b(2)*mu0);

% "multiply by x" operator
A = zeros([1,N-1]);
B = 0;
C = diag(ones(1,N-1));
D = zeros([N-1,1]);
X = cell2mat({A,B;C,D});

% coeff of the P_3 to P_N
disp(size(X));
for n = 2:N-1
   poly(:,n+1) = (X*poly(:,n)-a(n).*poly(:,n)-poly(:,n-1).*rtb(n-1))./rtb(n); 
end
end