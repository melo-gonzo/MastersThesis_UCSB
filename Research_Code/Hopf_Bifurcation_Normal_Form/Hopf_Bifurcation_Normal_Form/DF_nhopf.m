function [D, x_pos]=DF_nhopf(x_p,p)
if p.perturbed==0
    a=p.a(1);
else
    a=p.a(2);
end
b=p.b(1);
c=p.c(1);
d=p.d(1);
D=zeros(2,2);
x=x_p(1,1);
y=x_p(2,1);
x_pos=[x;y];
D(1,1)=a+2*x*(c*x-d*y)+c*(x^2+y^2);
D(1,2)=-b+2*y*(c*x-d*y)-d*(x^2+y^2);
D(2,1)=b+2*x*(d*x+c*y)+d*(x^2+y^2);
D(2,2)=a+2*y*(d*x+c*y)+c*(x^2+y^2);
end 