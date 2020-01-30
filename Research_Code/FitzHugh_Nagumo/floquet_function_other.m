function [kappa, V, scale_vector_idx]=floquet_function_other(p,s,numeric)
% kappa=-p.a(1)*2;
disp('floquet function other')

V=nan;
scale_vector_idx=nan;
if numeric==1
xg=s.xg;
tg=s.tg;
prc=s.prc(:,2:end);
% round(length(prc(:,1)/2))
prc_ic=prc(1,:);
prc_ic=[1;1];
% A=null(prc_ic);
% v1=A(:,1);
v1=[1; 0];

for k=1:length(s.xg(:,2))-1; if s.xg(k,2)<0 & s.xg(k+1,2)>0;sidx=k;continue;end;end
x0=[interp1(s.xg(sidx:sidx+1,2),s.xg(sidx:sidx+1,1),0);0];
% x0=xg(1,:);
% x0=[0;0];
w=0.1;
theta=(-1:0.01:1)*w;
% theta=[(-1:0.01:-0.5) 0.5:0.01:1]*w;
n1=x0(1)+v1(1)*theta;
n2=x0(2)+v1(2)*theta;
x_cross=nan(length(theta),2);
t_cross=nan(length(theta),1);
tspan=0:0.0001:tg(end)+0.1*tg(end);
% f0=[0 0];
for ii=1:length(theta)
    f0=[n1(ii);n2(ii)];
    crossLine=@(t,x) cross_line_other(t,x,prc_ic,f0);
    options=odeset('RelTol',1e-13,'AbsTol',1e-18*ones(length(x0),1),'Events',crossLine);
    [t,x,te,fe,~]=ode45(@(t,x) fhn_ode_equations(t,x,p),tspan,f0,options);
    x_cross(ii,:)=fe(end,:);
    te(end,:);
    t_cross(ii)=te(end);
    ii;
end
%% Brute force calculation
% delta_r=exp(-k*T);
dv=n1'-x0(1);
delta_r=abs(dv-(x_cross(:,1)-x0(1)));
fee=log(delta_r)./t_cross;
[dv delta_r fee]
% figure();plot(s.xg(:,1),s.xg(:,2));hold on;plot(n1,n2)
% plot(x(:,1),x(:,2))
% plot(fe(3,1),fe(3,2),'k.','markersize',20)
%%
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
U;
end 



