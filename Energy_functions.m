

% produced by Ozan Vahit Altınpınar (altinpinaro@itu.edu.tr) (2025)
% Matlab script Energy_functions.m
% Description: This script generates the graphs of the energy functions and probability density functions used in determining Similar Energy Regions.

clc;
clear all;

dmax = 10;
d=0:0.1:dmax;
N = length(d);
e = exp(1);
fun = @(x) x.*((1/(dmax*(exp(1)-1)))).*(exp(1-x/dmax));
T_H1 =  integral(fun,0,dmax)

fun = @(x) ((1/((e-1)*dmax))*(1-exp(1-x/dmax)) + (2*e-3)/((e-1)*dmax)).*x;
T_H2 =  integral(fun,0,dmax)

p_Nw = ((1/(dmax*(exp(1)-1)))).*(exp(1-d/dmax)); % probability density function for the narrow region
p_Mm(1:N) = 1/dmax; % probability density function for the medium region
p_We = ((1/(dmax*(exp(1)-1)))).*(1-(exp(1-(d)/dmax)))+(2*e-3)/((e-1)*dmax); % probability density function for the wide region

A_max = [max(p_Nw) max(p_Mm) max(p_We)];
p_max = max(A_max)+0.01;

t = 0:0.01:p_max;
M = length(t);
f_TH1(1:M) = T_H1;
f_TH2(1:M) = T_H2;


figure(1)
plot(d,p_We,'g','Linewidth',1)
hold on
plot(d,p_Mm,'b','Linewidth',1)
hold on
plot(d,p_Nw,'r','Linewidth',1)
hold on
plot(f_TH1,t,'y--','Linewidth',1)
hold on
plot(f_TH2,t,'m--','Linewidth',1)
hold on
xlabel("Measurement value (a) [m]  (a_{max} = 10 m) ")
ylabel("p_{proposed}(a)")
legend({'p(a|W_{e})','p(a|M_{m})','p(a|N_{w})','T_{H1}','T_{H2}'},'Location','southeast')

e_Nw = (1- (d/dmax).*exp(1-d/dmax)); % Energy function for the narrow region
e_Mm = (1-(d/dmax)); % Energy function for the medium region
e_We = (1 - ((d/((2*e-3)*dmax)).*(1-exp(1-d/dmax)) + d./dmax)); % Energy function for the wide region


figure(2)
plot(d,e_Mm,'b','Linewidth',1)
hold on
plot(d,e_Nw,'r','Linewidth',1)
hold on
plot(d,e_We,'g','Linewidth',1)

xlabel("Measurement value (a) [m]  (a_{max} = 10 m) ")
ylabel("e_{proposed}(a)")
legend({'e_{M_{m}}(a)','e_{N_{w}}(a)','e_{W_{e}}(a)'},'Location','southwest')
