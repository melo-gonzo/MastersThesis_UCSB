function [D, x_pos]=DF_thalamic(x_p,p)
if p.perturbed==0
    ib=p.ib(1);
else
    ib=p.ib(2);
end
cm=p.cm;
gl=p.gl;
el=p.el;
gna=p.gna;
ena=p.ena;
gk=p.gk;
ek=p.ek;
gt=p.gt;
et=p.et;
v=x_p(1);
h=x_p(2);
r=x_p(3);
x_pos=[v;h;r];
%dv
D(1,1)=cm.^(-1).*((-1).*gl+(-0.316406E0).*gk.*(1+(-1).*h).^4+(-1).*(1+ ...
    exp(1).^((1/7).*((-37)+(-1).*v))).^(-3).*gna.*h+(-1).*(1+exp(1).^( ...
  (-0.16129E0).*(60+v))).^(-2).*gt.*r+(-3/7).*exp(1).^((1/7).*((-37) ...
      +(-1).*v)).*(1+exp(1).^((1/7).*((-37)+(-1).*v))).^(-4).*gna.*h.*(( ...
      -1).*ena+v)+(-0.322581E0).*exp(1).^((-0.16129E0).*(60+v)).*(1+exp( ...
      1).^((-0.16129E0).*(60+v))).^(-3).*gt.*r.*((-1).*et+v));
  D(1,2)=cm.^(-1).*(0.126563E1.*gk.*(1+(-1).*h).^3.*((-1).*ek+v)+(-1).*(1+ ...
      exp(1).^((1/7).*((-37)+(-1).*v))).^(-3).*gna.*((-1).*ena+v));
  D(1,3)=(-1).*cm.^(-1).*(1+exp(1).^((-0.16129E0).*(60+v))).^(-2).*gt.*(( ...
      -1).*et+v);
  %dh
  D(2,1)=(-1/4).*exp(1).^((1/4).*(41+v)).*(1+exp(1).^((1/4).*(41+v))).^(-2) ...
      .*(0.128E0.*exp(1).^((1/18).*((-46)+(-1).*v))+4.*(1+exp(1).^((1/5) ...
      .*((-23)+(-1).*v))).^(-1))+((-0.711111E-2).*exp(1).^((1/18).*(( ...
      -46)+(-1).*v))+(4/5).*exp(1).^((1/5).*((-23)+(-1).*v)).*(1+exp(1) ...
      .^((1/5).*((-23)+(-1).*v))).^(-2)).*((1+exp(1).^((1/4).*(41+v))) ...
      .^(-1)+(-1).*h);
  D(2,2)=(-0.128E0).*exp(1).^((1/18).*((-46)+(-1).*v))+(-4).*(1+exp(1).^(( ...
      1/5).*((-23)+(-1).*v))).^(-1);
  D(2,3)=0;
  %dr
  D(3,1)=(-1/4).*exp(1).^((1/4).*(84+v)).*(28+exp(1).^((-0.952381E-1).*(25+ ...
      v))).^(-1).*(1+exp(1).^((1/4).*(84+v))).^(-2)+0.952381E-1.*exp(1) ...
      .^((-0.952381E-1).*(25+v)).*(28+exp(1).^((-0.952381E-1).*(25+v))) ...
      .^(-2).*((1+exp(1).^((1/4).*(84+v))).^(-1)+(-1).*r);
  D(3,2)=0;
  D(3,3)=(-1).*(28+exp(1).^((-0.952381E-1).*(25+v))).^(-1);
end





