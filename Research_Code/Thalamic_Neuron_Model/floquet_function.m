function [kappa, V, scale_vector_idx]=floquet_function(p,s,numeric)
disp('Calculating Floquet Exponents')
% kappa=-p.a(1)*2;
V=nan;
scale_vector_idx=nan;
if numeric==1
xg=s.xg;
tg=s.tg;
prc=s.prc(:,2:end);
% round(length(prc(:,1)/2))
prc_ic=prc(1,:);
A=null(prc_ic);
v1=A(:,1);
v2=A(:,2);
x0=xg(1,:);
w=0.00001;
theta=0:0.01*2*pi:2*pi;
% theta=[(-1:0.1:-0.5) 0.5:0.1:1]*w;
n1=x0(1)+w*(v1(1).*cos(theta)+v2(1).*sin(theta));
n2=x0(2)+w*(v1(2).*cos(theta)+v2(2).*sin(theta));
n3=x0(3)+w*(v1(3).*cos(theta)+v2(3).*sin(theta));
x_cross=nan(length(theta),3);
t_cross=nan(length(theta),1);
tspan=0:p.dt:tg(end)+0.1*tg(end);
% f0=[0 0];
for ii=1:length(theta)
    f0=[n1(ii);n2(ii);n3(ii)];
    crossplane=@(t,x) cross_plane(t,x,prc_ic,f0);
    options=odeset('RelTol',3e-14,'AbsTol',1e-17*ones(length(x0),1),'Events',crossplane);
    [t,x,te,fe,~]=ode45(@(t,x) thalamic_ode_equations(t,x,p),tspan,f0,options);
    x_cross(ii,:)=fe(end,:);
    te(end,:);
    t_cross(ii)=te(end);
    ii;
end
x0=[n1;n2;n3];
% DP=(x_cross-xg(1,:))'/(x0-xg(1,:)');
DP=(x_cross)'/(x0);
% DP=x_cross'/x0; 
[V,U]=eig(DP);
%Scaling of eigenvectors to have unit length
for ii=1:3
    if round(sum(U(:,ii)),3)~=1
        V(:,ii)=V(:,ii)-prc_ic*V(:,ii)*prc_ic'/(norm(prc_ic)^2);
    end
end
k=log(U)/tg(end); %Floquet exponents
if ~isreal(k)
    disp('imaginary')
end
kappa=sort(diag(k));
% disp('check scale_vector_idx floquet')
scale_vector_idx=find(diag(k)==min(diag(k)));
U;
end 



