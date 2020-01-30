% close all
% clear all
% clc
% set_default_plot
t=0:0.001:20;
c=-1;
d=1;
b=1;
a=0.1;

q=sqrt(-a/c);


theta=b*t-(d/(2*c))*log(a+c*q^2-c*q^2*exp(2*a*t))+(d/(2*c))*log(a);
% plot(theta)


x=cos((b-d*a/c)*t).*(sqrt(a)*exp(a*t)*q./sqrt((a+c*q^2-c*exp(2*a*t)*q^2)));
y=sin((b-d*a/c)*t).*(sqrt(a)*exp(a*t)*q./sqrt((a+c*q^2-c*exp(2*a*t)*q^2)));

% plot(x);axis tight
% figure();
% sx=exp(a*t).*(2*a*t+exp(2*a*t)+1).*cos(t)./(2*sqrt(a).*(exp(2*a*t)+1).^(3/2));
% plot(sx);axis tight

% sx=(exp(a*t).*q.*(-c*(-2*a^2*t+c*q^2*(-1+exp(2*a*t)-2*a*t)).*cos((b-a*d/c)*t)+...
%     2*a*d.*(a-c*(-1+exp(2*a*t))*q^2).*t.*sin((b-a*d/c)*t)))./...
%     (2*sqrt(a)*c.*(a-c*(-1+exp(2*a*t))*q^2).^(3/2));
% 
% sy=(exp(a*t).*q.*(-c*(-2*a^2*t+c*q^2*(-1+exp(2*a*t)-2*a*t)).*sin((b-a*d/c)*t)-...
%     2*a*d.*(a-c*(-1+exp(2*a*t))*q^2).*t.*cos((b-a*d/c)*t)))./...
%     (2*sqrt(a)*c.*(a-c*(-1+exp(2*a*t))*q^2).^(3/2));
shit=(2*b*c*t+d*log(a)-d*log(a-c*(-1+exp(2*a*t))*q^2))/(2*c);


sx=(exp(a*t).*q.*((2*a^2.*t-c*q^2.*(-1+exp(2*a*t)-2*a*t)).*cos(shit)-d*q^2.*(1+exp(2*a*t).*(-1+2*a.*t)).*sin(shit)))./(2*sqrt(a)*(a-c*(-1+exp(2*a*t)).*q^2).^(3/2));

sy=(exp(a*t).*q.*((2*a^2.*t-c*q^2.*(-1+exp(2*a*t)-2*a*t)).*sin(shit)+d*q^2.*(1+exp(2*a*t).*(-1+2*a.*t)).*cos(shit)))./(2*sqrt(a)*(a-c*(-1+exp(2*a*t)).*q^2).^(3/2));



% 
% plot(sx);axis tight
% hold on;plot(sy)