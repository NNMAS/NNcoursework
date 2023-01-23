%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Problem 1 OU-LIF 
%  Jan 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; 

%parameters
t_start = 0;          %simulation start time
t_end = 1000;          %simuation end time
dt = 0.01;            %time step
tau = 0.1;            %relaxation time
sigma=1;
mu=0.2;
x0 = 1;               
m1 = 0;               
y0 = 0;               
start_dist = -2.0;    
end_dist = 2.0;       

%time
T = t_start:dt:t_end;

% compute x and y
i = 1;
x(1) = x0; 
y(1) = y0;
for t=t_start+dt:dt:t_end
   i = i + 1; 
   r1 = randn;
   r2 = randn;
   x(i) = x(i-1)*exp(-dt/tau) +mu*dt+ sqrt((tau*0.5)*(1-(exp(-dt/tau))^2))*sigma*r1;
   y(i) = y(i-1) + x(i-1)*tau*(1-exp(-dt/tau))+sigma*sqrt((tau^3*(dt/tau-2*(1-exp(-dt/tau))+0.5*(1-exp(-2*dt/tau))))-((0.5*sigma^2*tau^2)*(1-exp(-dt/tau))^2)^2/((sigma^2*tau/2)*(1-exp(-2*dt/tau))))*r2+((0.5*sigma^2*tau^2)*(1-exp(-dt/tau))^2)/(sigma*sqrt((tau/2)*(1-(exp(-dt/tau))^2)))*r1;
end

% pdf for OU process
k = 0; j = start_dist:dt:end_dist;
for l=start_dist:dt:end_dist
    k = k + 1;
    p(k) = (sqrt((1/tau)/pi)/sigma)*exp(-(1/tau)*(l-m1)^2/(sigma)^2); 
end

mean(x)
var(x)

figure;
subplot(3,1,1)
plot(T,x,'k-')
xlabel('time')
ylabel('x')
subplot(3,1,2)
plot(T,y,'k-')
xlabel('time')
ylabel('y')
subplot(3,1,3)
hold on
histogram(x,60)
plot(j,p,'r-')
xlabel('x')
ylabel('probability')
hold off