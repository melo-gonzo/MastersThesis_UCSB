function PSI=a_approach(params,r)
% input r=nan to calculate psi from params, or give specific r to calc psi
% at, or a vector to give full psi
a=params.a(1);
b=params.b(1);
c=params.c(1);
d=params.d(1);
if isnan(r)
    r=sqrt(-a/c);
end

c1=-sqrt((-a*c*(c^2+d^2)))/2;
PSI=c1*(a/r.^2 +c);
