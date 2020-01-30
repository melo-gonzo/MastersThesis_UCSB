function xout=thalamic_ode_equations(~,x,p)
cm=p.cm;
gl=p.gl;
el=p.el;
gna=p.gna;
ena=p.ena;
gk=p.gk;
ek=p.ek;
gt=p.gt;
et=p.et;

ib=p.ib;
hinf=1/(1+exp((x(1)+41)/4));
rinf=1/(1+exp((x(1)+84)/4));
alphah=0.128*exp(-(x(1)+46)/18);
betah=4/(1+exp(-(x(1)+23)/5));
tauh=1/(alphah+betah);
taur=28+exp(-(x(1)+25)/10.5);
minf=1/(1+exp(-(x(1)+37)/7));
pinf=1/(1+exp(-(x(1)+60)/6.2));
iL=gl*(x(1)-el);
ina=gna*(minf^3)*x(2)*(x(1)-ena);
ik=gk*((0.75*(1-x(2)))^4)*(x(1)-ek);
it=gt*(pinf^2)*x(3)*(x(1)-et);

xout=[(-iL-ina-ik-it+ib)/cm;
    (hinf-x(2))/tauh;
    (rinf-x(3))/taur];

