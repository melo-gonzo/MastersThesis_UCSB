function theta = a_phase(p,r,phi)
%asymptotic phase
if isnan(r)
    r=sqrt(-p.a(1)/p.c(1));
end 
if isnan(phi)
    phi=0;
end 
c2=(p.d(1)/(2*p.c(1)))*log(-p.a(1)/p.c(1));
theta=phi-(p.d(1)/p.c(1))*log(r)+c2;
% theta=phi+(p.d(1)/(2*p.c(1)))*log(p.a(1))-(p.d(1)/(2*p.c(1)))*log(-p.c(1)*r.^2);
