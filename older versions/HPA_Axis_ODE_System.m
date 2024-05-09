%Initialize the variables
t0 = 0;
tf = 48;
x = [120, 1.2, 0.0065, 0];

%Define the range of time and the initial conditions
tspan = [t0, tf];

%Solve the ODE
[t,x] = ode45(@HPA_Axis_ODE_System, tspan, x0);

%Plot the results
plot(t,x(:,1),'b-',t,x(:,2),'r-',t,x(:,3),'g-',t,x(:,4),'m-');
title('HPA Axis ODE System');
xlabel('Time (hr)');
ylabel('Concentration (uM)');
legend('CRH','ACTH','GR','Cortisol');

function dAlldt = HPA_Axis_ODE_System(t,x)

%intialize constants
kcd = 1;        %CRH degradation constant

kad = 10;       %ACTH degradation constant

krd = 0.9;      %GR degradation constant 

kcr = 0.05;     %GR production constant

k = 0.001;      %equilibrium binding affinity

ki1 = 0.1;      %inhibition constant 1

ki2 = 0.1;      %inhibtion constant 2

Kc = 1.0;       

Ka = 1.0;       %ACTH production constant

Kr = 1.0;      

Kod = 1.0;      %Cortisol degradation constant

Ko = 1.0;       %Cortisol production constant

C = x(1);       %initial value of CRH
A = x(2);       %intial value of ACTH
R = x(3);       %initial value of GR
O = x(4);       %intial value of Cortisol
or = x(3);      %initial value of dimerized glucocorticoid receptor

f = 0.02+0.01*sin(2*pi*t/12);

%define terms of differential equations
c = ((Kod*C)/Kc);
a = ((Kod^2)*A)/Kc*Ka;
o = ((Kod^3)*O)/(Kc*Ka*Ko);
r = ((Kod*R)/Kr);


%define differential equations
dcdt = (1+f/(1+(o/ki1)))- kcd*c;
dadt = (c/(1+(or/ki2)))-kad*a;
drdt = (or^2)/(k+or^2)+kcr-krd*r;
dodt = a - o;

%define function as system of differential equations
dAlldt = [dcdt; dadt; drdt; dodt];