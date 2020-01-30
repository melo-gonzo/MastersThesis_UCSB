function [kappa, V, scale_vector_idx]=floquet_function(p,s,numeric)
kappa=-p.a(1)*2;
V=nan;
scale_vector_idx=nan;
if numeric==1
xg=s.xg;
tg=s.tg;
prc=s.prc(:,2:end);

prc_ic=prc(1,:);
A=null(prc_ic);
v1=A(:,1);
x0=xg(1,:);
w=0.0001;
theta=(-1:0.1:1)*w;
n1=x0(1)+v1(1)*theta;
n2=x0(2)+v1(2)*theta;
x_cross=nan(length(theta),2);
t_cross=nan(length(theta),1);
tspan=0:p.dt:tg(end)+0.1*tg(end);
% f0=[0 0];
for ii=1:length(theta)
    f0=[n1(ii);n2(ii)];
    crossLine=@(t,x) cross_line(t,x,prc_ic,f0);
    options=odeset('RelTol',3e-14,'AbsTol',1e-17*ones(length(x0),1),'Events',crossLine);
    [t,x,te,fe,~]=ode45(@(t,x) hopf_ode_equations(t,x,p),tspan,f0,options);
    x_cross(ii,:)=fe(end,:);
    te(end,:);
    t_cross(ii)=te(end);
    ii;
end
x0=[n1;n2];
DP=x_cross'/x0;

[V,U]=eig(DP);
%Scaling of eigenvectors to have unit length
for ii=1:2
    if round(sum(U(:,ii)),3)~=1
        V(:,ii)=V(:,ii)-prc_ic*V(:,ii)*prc_ic'/(norm(prc_ic)^2);
    end
end
k=log(U)/tg(end); %Floquet exponents
if ~isreal(k)
    disp('imaginary')
end
kappa=sort(diag(k));
scale_vector_idx=find(diag(k)==min(diag(k)));
V(:,scale_vector_idx)=-V(:,scale_vector_idx);
U;
end 



% DP=(x_cross-xg(1,:))'/(x0-xg(1,:)');
% DP=(x0'\x_cross)'