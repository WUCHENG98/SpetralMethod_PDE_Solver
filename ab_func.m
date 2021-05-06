% compute the recurrence coefficients alpha and beta in Guass Stieltjes procedure


function ab_func(func,name,n,xmin,xmax,nint,npts)
% This function uses Guass Stieltjes procedure to calculate recurrence coefficients alpha and beta ...
% for weight function w(x)=func(x). The integral is on [xmin,xmax]; 
% which is divdided into subintervals with number of points per interval  
% also mu_0 would be saved in the third column

% input:
%    func:   the function handle of the weight function
%    name:   the name of the weight function used for save the .dat file
%    m:      the number of alpha and beta needed
%    xmin:   left end point of the integral interval
%    xmax:   right end point of the integral interval
%    nint:   the number of interval
%    npts:   the number of points per interval

% Output:
%    no returned value; alpha's and beta's are saved in the .dat file 
%    (also mu_0 would be saved in the third column)

format long e

% find quadrature pts and wts and weight functions
[pts,wts] = multidomquad(nint,npts,xmin,xmax); 
wtfcn=func(pts);


% construct polynomial Q_0(x)
ntot=nint*npts; 
Q0=ones(ntot,1);

% first two integrals
s1=sum(wts.*wtfcn); 
s2=sum(wts.*(pts.*wtfcn));

% find norm and alpha_0 and "beta_0"
h(1)=s1; 
a(1)=s2/s1; 
b(1)=0;

% store alpha and beta (and mu0 = s1)
fileName= ['ab_',name,'.dat'];
myfile = fopen(fileName, 'wt');
fprintf(myfile,'%20.12f %20.12f %20.12f\n',a(1),b(1),s1);

% find polynomial Q_1(x)
Q1=pts-a(1);
% next two integrals 
s1=sum(wts.*(wtfcn.*(Q1.^2))); 
s2=sum(wts.*(pts.*(wtfcn.*(Q1.^2))));

% norm and alpha_1 and beta_1
a(2)=s2/s1; 
h(2)=s1; 
b(2)=h(2)/h(1);

fprintf(myfile,'%20.12f %20.12f\n',a(2),b(2));

for k=3:n
    % find the next polynomial
    Q2=(pts-a(k-1)).*Q1-b(k-1)*Q0;
    % find the next integral
    s1=sum(wts.*(wtfcn.*(Q2.^2)));
    s2=sum(wts.*(pts.*(wtfcn.*(Q2.^2))));
    % find the next alpha, norm and beta
    a(k)=s2/s1; h(k)=s1; b(k)=h(k)/h(k-1);
    % store
    fprintf(myfile,'%20.12f %20.12f\n',a(k),b(k));
    % change index of polys
    Q0=Q1; Q1=Q2;
end
end


function [pts,wts] = multidomquad(nint,npts,xmin,xmax)
% Generate the Fejer quadraturesa and weights on [xmin,xmax] which is divdided into subintervals 
% with number of points per interval

% Input:
%    nint:   the number of interval
%    npts:   the number of points per interval
%    xmin:   left end point of the integral interval
%    xmax:   right end point of the integral interval

% Output:
%    pts and wts as required
   
format long e        
ntot = nint*npts; % total points
dx = (xmax-xmin)/ntot;

pts = [];
wts = [];

for i=1:nint
    % find subinterval [a,b]
    a = xmin+(i-1)*npts*dx;  
    b = a+npts*dx; 

    [subpts,subwts]=fejer(a,b,npts); 
    pts = [pts subpts'];
    wts = [wts subwts'];
end
pts = pts';
wts = wts';
end




function [subpts,subwts]=fejer(a,b,npts) 
% Generates the n-point quadrature in [a,b]

% Input:
%   a:      left end of the interval
%   b:      right end of the interval
%   npts:   the number of points in the interval

% Output:
%   subpts and subwts of the subinterval

format long e

% Fejer quadrature pts x_i = cos[(2i-1)pi/(2N)]
i = npts:-1:1; 
theta = (2*i-1)*pi./(2*npts); 
x = cos(theta');  

% Fejer quadrature wts w_i = (2/N)[1-sum_{j=1}^{N/2} cos(2jx_i)/(4j^2-1)]
j=1:floor(npts/2); 
for i_loop=npts:-1:1 
  s=sum(cos(2*j*theta(i_loop))./(4*(j.^2)-1)); 
  w(i_loop)=2*(1-2*s)/npts; 
end

% mapping [-1,1] onto [a,b]. 
r1 = (b-a)/2.0; 
r2 = r1+a; 
subpts = r1*x+r2; 
subwts = (r1*w)'; 
end