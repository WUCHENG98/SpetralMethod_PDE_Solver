function [pts,wts] = GSprocedure(name,N)
% find the quadrature points and weigths with GS procedure

% Input: 
%    name: name of the .dat file which stores the alphas and betas
%    N:    the number of basis / quadrature pts required

% Output:
%   pts,wts: quadrature points and weigths

format long e

% load alpha and beta (and mu0)
fileName= ['ab_',name,'.dat'];
data = importdata(fileName);
a = data(:,1);
b = data(:,2);
mu0 = data(1,3);
rtb = sqrt(b);
rtb(1)=[];

% construct Jacobi matrix and find eigenvectors/values
J=diag(rtb(1:N-1),-1)+diag(a(1:N))+diag(rtb(1:N-1),1);
[v,e] = eig(J);
pts = diag(e);
wts = mu0.*v(1,:)'.^2;
end