%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Problem 4 E-I neural network
%  Jan 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% parameters
tau_I=1:0.05:3;
h_E=0.6;
tau_E=0.3;
M_EE=0.5;
M_EI=0.5;
M_IE=0.1;
M_II=0.2;
h_I=0.4;
tspan = [0 10];
v0 = [1 1];
vE=[];
vI=[];
for i=1:length(tau_I)
    [t,v] = ode45(@(t,v) EIdyna(t,v,tau_E,tau_I(i),M_EE,M_EI,M_IE,M_II,h_E,h_I),tspan,v0);
    plot(t,v(:,1),'-o',t,v(:,2),'-*','LineWidth',1.5)
    hold on;
    vE=[vI;v(end,1)];
    vI=[vI;v(end,2)];
end

% Figure 1
xlabel('Time t', 'Interpreter','latex','FontSize',14);
ylabel('Solution', 'Interpreter','latex','FontSize',14);
title('Solution of E-I neural network with $\tau_{I}\in(1,3)$', 'Interpreter','latex','FontSize',16);

% Figure 2
figure();
plot(h_E,vE,'-o',h_E,vI,'-*','LineWidth',1.5)


function dydt = EIdyna(t,v,tau_E,tau_I,M_EE,M_EI,M_IE,M_II,h_E,h_I)
 dydt = zeros(2,1);
 dydt(1) = (-v(1)+max(M_EE*v(1)-M_EI*v(2)+h_E,0))/tau_E;
 dydt(2) = (-v(2)+max(M_IE*v(1)-M_II*v(2)+h_I,0))/tau_I;
end