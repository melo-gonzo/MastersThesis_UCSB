function xout=hopf_ode_equations(~,x,p)
xout=[p.a(1)*x(1)-p.b(1)*x(2)+(p.c(1)*x(1)-p.d(1)*x(2))*(x(1)^2+x(2)^2);
    p.b(1)*x(1)+p.a(1)*x(2)+(p.d(1)*x(1)+p.c(1)*x(2))*(x(1)^2+x(2)^2)];