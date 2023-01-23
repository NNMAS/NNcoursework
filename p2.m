%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Problem 2 Hodgkin-Huxley Model
%  Jan 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

%parameters
vRest = -56;
v = -20:0.01:60;
u = vRest - v;
alpha_m = -0.08 * (56+u) ./ (exp(-(56+u)/6.8) - 1);
beta_m = 0.8 * exp(-(56+u)/18);

alpha_h =  0.006 * exp(-(41+u)/14.7);
beta_h =  1.3 ./ (exp(-(41+u)/7.6) + 1);

alpha_n =  0.0088 * (-(40+u)) ./ (exp(-(40+u)/7) - 1);
beta_n =  0.037 * exp(-(40+u)/40);

% time constants
tau_n = 1 ./ (alpha_n + beta_n);
tau_m = 1 ./ (alpha_m + beta_m);
tau_h = 1 ./ (alpha_h + beta_h);

figure;
hold on; grid minor;
plot(v, tau_h, 'LineWidth', 2)
plot(v, tau_m, 'LineWidth', 2)
plot(v, tau_n, 'LineWidth', 2)
title('time constants', 'Interpreter','latex')
legend('\tau_h', '\tau_m', '\tau_n')

% steady state values
n_inf = alpha_n ./ (alpha_n + beta_n);
m_inf = alpha_m ./ (alpha_m + beta_m);
h_inf = alpha_h ./ (alpha_h + beta_h);

figure;
hold on; grid minor;
plot(v, h_inf, 'LineWidth', 2)
plot(v, m_inf, 'LineWidth', 2)
plot(v, n_inf, 'LineWidth', 2)
title('steady state values', 'Interpreter','latex')
legend('$$h_{\infty}$$', '$$m_{\infty}$$', '$$n_{\infty}$$', 'Interpreter', 'latex')


%% Current Test
% very low input current

dt = 0.05; % Simulation time step
Duration = 4000; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms

I = 2*ones(1, T);
[v, n, m, h] = simulate_HH(dt, I, T);

figure;
subplot(1, 3, 1)
plot(v, n); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('n', 'Interpreter','latex')
subplot(1, 3, 2)
plot(v, m); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('m', 'Interpreter','latex')
subplot(1, 3, 3)
plot(v, h); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('h', 'Interpreter','latex')
sgtitle('$$I = 2\mu A$$', 'Interpreter', 'latex')

% normal input current
I = 8*ones(1, T);
[v, n, m, h] = simulate_HH(dt, I, T);

figure;
subplot(1, 3, 1)
plot(v, n); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('n', 'Interpreter','latex')
subplot(1, 3, 2)
plot(v, m); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('m', 'Interpreter','latex')
subplot(1, 3, 3)
plot(v, h); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('h', 'Interpreter','latex')
sgtitle('$$I = 8\mu A$$', 'Interpreter', 'latex')

% very high input current
I = 200*ones(1, T);
[v, n, m, h] = simulate_HH(dt, I, T);

figure;
subplot(1, 3, 1)
plot(v, n); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('n', 'Interpreter','latex')
subplot(1, 3, 2)
plot(v, m); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('m', 'Interpreter','latex')
subplot(1, 3, 3)
plot(v, h); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('h', 'Interpreter','latex')
sgtitle('$$I = 200\mu A$$', 'Interpreter', 'latex')

%%
function [v, n, m, h] = simulate_HH(dt, I, T)
% constants
Cm = 1.9; % Membrane capacitance in micro Farads
gNa = 50; % in Siemens, maximum conductivity of Na+ Channel
gK = 22; % in Siemens, maximum conductivity of K+ Channel
gl = 0.4; % in Siemens, conductivity of leak Channel
ENa = 50; % in mv, Na+ nernst potential
EK = -70; % in mv, K+ nernst potential
El = -81; % in mv, nernst potential for leak channel
vRest = -56; % in mv, resting potential


alpha_m = @(u) -0.08 * (56+u)./(exp(-(56+u)/6.8) - 1);
beta_m = @(u) 0.8 * exp(-(56+u)/18);

alpha_h = @(u) 0.006 * exp(-(41+u)/14.7);
beta_h = @(u) 1.3 ./ (exp(-(41+u)/7.6) + 1);

alpha_n = @(u) 0.0088 * (-(40+u)) ./ (exp(-(40+u)/7) - 1);
beta_n = @(u) 0.037 * exp(-(40+u)/40);


% initial values
v = zeros(1, T); % output voltage
v(1) = vRest;
n = zeros(1, T); % probability of K+  activation gate being open
n(1) = alpha_n(vRest - v(1)) / (alpha_n(vRest - v(1)) + beta_n(vRest - v(1)));
m = zeros(1, T); % probability of Na+ activation gate being open
m(1) = alpha_m(vRest - v(1)) / (alpha_m(vRest - v(1)) + beta_m(vRest - v(1)));
h = zeros(1, T); % probability of NA+ inactivation gate being open
h(1) = alpha_h(vRest - v(1)) / (alpha_h(vRest - v(1)) + beta_h(vRest - v(1)));
% 2nd order Runge-Kutta
k1 = zeros(4, 1);
k2 = zeros(4, 1);
for i = 1:T-1
    k1(1) = dt * (-1)/Cm * (gl*(v(i)-El) + gK*n(i)^4*(v(i)-EK) + gNa*m(i)^3*h(i)*(v(i)-ENa)-I(i));
    tau_n = 1 / (alpha_n(vRest - v(i)) + beta_n(vRest - v(i)));
    n_inf = alpha_n(vRest - v(i)) / (alpha_n(vRest - v(i)) + beta_n(vRest - v(i)));
    k1(2) = dt * 1/tau_n * (-n(i) + n_inf);
    tau_m = 1 / (alpha_m(vRest - v(i)) + beta_m(vRest - v(i)));
    m_inf = alpha_m(vRest - v(i)) / (alpha_m(vRest - v(i)) + beta_m(vRest - v(i)));
    k1(3) = dt * 1/tau_m * (-m(i) + m_inf);
    tau_h = 1 / (alpha_h(vRest - v(i)) + beta_h(vRest - v(i)));
    h_inf = alpha_h(vRest - v(i)) / (alpha_h(vRest - v(i)) + beta_h(vRest - v(i)));
    k1(4) = dt * 1/tau_h * (-h(i) + h_inf);

    k2(1) = dt * (-1)/Cm * (gl*((v(i)+k1(1))-El) + gK*(n(i)+k1(2))^4*((v(i)+k1(1))-EK) + gNa*(m(i)+k1(3))^3*(h(i)+k1(4))*((v(i)+k1(1))-ENa)-I(i));
    tau_n = 1 / (alpha_n(vRest - (v(i)+k1(1))) + beta_n(vRest - (v(i)+k1(1))));
    n_inf = alpha_n(vRest - (v(i)+k1(1))) / (alpha_n(vRest - (v(i)+k1(1))) + beta_n(vRest - (v(i)+k1(1))));
    k2(2) = dt * 1/tau_n * (-(n(i)+k1(2)) + n_inf);
    tau_m = 1 / (alpha_m(vRest - (v(i)+k1(1))) + beta_m(vRest - (v(i)+k1(1))));
    m_inf = alpha_m(vRest - (v(i)+k1(1))) / (alpha_m(vRest - (v(i)+k1(1))) + beta_m(vRest - (v(i)+k1(1))));
    k2(3) = dt * 1/tau_m * (-(m(i)+k1(3)) + m_inf);
    tau_h = 1 / (alpha_h(vRest - (v(i)+k1(1))) + beta_h(vRest - (v(i)+k1(1))));
    h_inf = alpha_h(vRest - (v(i)+k1(1))) / (alpha_h(vRest - (v(i)+k1(1))) + beta_h(vRest - (v(i)+k1(1))));
    k2(4) = dt * 1/tau_h * (-(h(i)+k1(4)) + h_inf);

    v(i+1) = v(i) + (k1(1)+k2(1))/2;
    n(i+1) = n(i) + (k1(2)+k2(2))/2;
    m(i+1) = m(i) + (k1(3)+k2(3))/2;
    h(i+1) = h(i) + (k1(4)+k2(4))/2;
end

end







