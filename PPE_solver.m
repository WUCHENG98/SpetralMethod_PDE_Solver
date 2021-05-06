function [eigval,eigfunc]=PPE_solver(name,N)
% solve the FPE with L operator defined inside based on the basis set with
% given name of data file, the number of basis polynomials and the weight
% function.

format long e

% define B(x)
s=1;
g = sqrt(0.05);
h = @(gm)(3*sqrt(pi)/4)/gm^3;
G1 = @(x,gm) h(gm)*(erf(x.*gm) - 2./sqrt(pi).*x.*gm.*exp(-x.^2.*gm^2));
B = @(x) G1(s.*x,g)./(s^2.*x.^3);

% pts and wts
[pts,wts] = GSprocedure(name,N);

% generate the polymatrix
poly = polygen(name,N);

% find the derivative operator
D = polydif(pts,wts,poly);

% define the L operator
L = D'*diag(B(pts))*D;
%*diag(B(pts))

% find the eigenvalues and eigenfunctions and sort them
[v,e]=eig(L);
[eigval, ind] = sort(diag(e));
eigfunc = v(:,ind);

disp('eigenvalues')
for i=1:N
fprintf('%2i %20.16e\n', i-1, eigval(i))
end

for i=1:N
    indexedeigf2(i)=eigfunc(i,1); %MATLAB has no zero index
    indexedeigf4(i)=eigfunc(i,2);
    indexedeigf6(i)=eigfunc(i,3);
    indexedeigf8(i)=eigfunc(i,4);
    indexedeigf10(i)=eigfunc(i,5);
    indexedeigf12(i)=eigfunc(i,6);
    indexedeigf12(i)=eigfunc(i,7);
end



figure(1);
subplot(2,3,1);
%plot the numeig'th eigenfunction
plot(pts,indexedeigf2,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([0 15 -1 1]) %adjust axis here, try xmin and 15
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.5:1],'linewidth',1.6)
set(gca,'Xtick',[0:5:15],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_0(x)$','Interpreter','Latex','fontsize',20)


subplot(2,3,2);
plot(pts,indexedeigf4,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([0 15 -1 1]) %adjust axis here, try xmin and 15
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.5:1],'linewidth',1.6)
set(gca,'Xtick',[0:5:15],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_1(x)$','Interpreter','Latex','fontsize',20)

subplot(2,3,3);
plot(pts,indexedeigf6,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([0 15 -1 1]) %adjust axis here, try xmin and 15
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.5:1],'linewidth',1.6)
set(gca,'Xtick',[0:5:15],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_2(x)$','Interpreter','Latex','fontsize',20)

subplot(2,3,4);
plot(pts,indexedeigf8,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([0 15 -1 1]) %adjust axis here, try xmin and 15
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.5:1],'linewidth',1.6)
set(gca,'Xtick',[0:5:15],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_3(x)$','Interpreter','Latex','fontsize',20)

subplot(2,3,5);
plot(pts,indexedeigf10,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([0 15 -1 1]) %adjust axis here, try xmin and 15
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.5:1],'linewidth',1.6)
set(gca,'Xtick',[0:5:15],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_{4}(x)$','Interpreter','Latex','fontsize',20)

subplot(2,3,6);
plot(pts,indexedeigf12,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([0 15 -1 1]) %adjust axis here, try xmin and 15
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.5:1],'linewidth',1.6)
set(gca,'Xtick',[0:5:15],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_{5}(x)$','Interpreter','Latex','fontsize',20)


end

