% produced by Ozan Vahit Altınpınar (altinpinaro@itu.edu.tr) (2025)
% Matlab script Energy_functions.m
% Description: This script generates the graphs of the energy functions and probability density functions used for determining Similar Energy Regions.

clc;
clear all;

amax = 10;
a=0:0.1:amax;
N = length(a);
e = exp(1);
fun = @(x) x.*((1/(amax*(e-1)))).*(exp(1-x/amax));
T_H1 =  integral(fun,0,amax)

fun = @(x) ((1/((e-1)*amax))*(1-exp(1-x/amax)) + (2*e-3)/((e-1)*amax)).*x;
T_H2 =  integral(fun,0,amax)

p_Nw = ((1/(amax*(e-1)))).*(exp(1-a/amax)); % probability density function for the narrow region
p_Mm(1:N) = 1/amax; % probability density function for the medium region
p_We = ((1/(amax*(e-1)))).*(1-(exp(1-(a)/amax)))+(2*e-3)/((e-1)*amax); % probability density function for the wide region

A_max = [max(p_Nw) max(p_Mm) max(p_We)];
p_max = max(A_max)+0.01;

t = 0:0.01:p_max;
M = length(t);
f_TH1(1:M) = T_H1;
f_TH2(1:M) = T_H2;


figure(1)
plot(a,p_We,'g','Linewidth',1)
hold on
plot(a,p_Mm,'b','Linewidth',1)
hold on
plot(a,p_Nw,'r','Linewidth',1)
hold on
plot(f_TH1,t,'y--','Linewidth',1)
hold on
plot(f_TH2,t,'m--','Linewidth',1)
hold on
xlabel("Measurement value (a) [m]  (a_{max} = 10 m) ")
ylabel("p_{proposed}(a)")
title("Probability Density Functions")
legend({'p(a|W_{e})','p(a|M_{m})','p(a|N_{w})','T_{H1}','T_{H2}'},'Location','southeast')

e_Nw = (1- (a/amax).*exp(1-a/amax)); % Energy function for the narrow region
e_Mm = (1-(a/amax)); % Energy function for the medium region
e_We = (1 - ((a/((2*e-3)*amax)).*(1-exp(1-a/amax)) + a./amax)); % Energy function for the wide region


figure(2)
plot(a,e_Mm,'b','Linewidth',1)
hold on
plot(a,e_Nw,'r','Linewidth',1)
hold on
plot(a,e_We,'g','Linewidth',1)

xlabel("Measurement value (a) [m]  (a_{max} = 10 m) ")
ylabel("e_{proposed}(a)")
title("Energy Functions")
legend({'e_{M_{m}}(a)','e_{N_{w}}(a)','e_{W_{e}}(a)'},'Location','southwest')
