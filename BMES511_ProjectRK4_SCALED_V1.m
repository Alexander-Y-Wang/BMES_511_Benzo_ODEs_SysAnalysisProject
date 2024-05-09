clear all
%% Initial 
C = .05; %corticotropin releasing hormone
A = .5; %adrenocorticotropic hormone
O = 0.25; %cortisol
R = .06;%glucocorticoid receptors on adrenal

%% Constants
%CRK (C): Corticotropin releasing hormone
Kc = 1; %Production of Corticotropin releasing hormone
Kcd = 1; %Constant degradation of Corticotropin releasing hormone
Kor = 0.05; %pituitary GR production

% a (ACTH): adrenocorticotropic hormone
Ka = 1; %Production of adrenocorticotropic hormone
Kad = 10; %Constant degradation of adrenocorticotropic hormone

% o (O): Cortisol
Ko = 1; %Production of Cortisol
Kod = 1; %Constant degradation of Cortisol

% r (R): glucocorticoid receptor (IN PITIUTARY) 
Kr = 1; %Production of glucocorticoid receptor
Krd = .8;  %Constant degradation of glucocorticoid receptor !!!!!!!!!!!!

% i : inhibition constants
Ki1 = 0.1; %inhibition 1
Ki2 = 0.1; %inhibition 2

K = 0.001; %equilibrium binding affinity

T = [0]; %initial time
dt = 0.0001; %time step
T_end = 10; %time end


%SCALED!!!!!!!!!!
c = (Kod*C)/Kc;
a = ((Kod^2)*A)/(Kc*Ka);
o = ((Kod^3)*O)/(Kc*Ka*Ko);
r = ((Kod*R)/Kr);
kcd = Kcd/Kod;
kad = Kad/Kod;
krd = Krd/Kod;

for i = 1:(T_end/dt)
    if i <= 10000
        F = 100; %This stress value was pulled from Srisam (PTSD paper)
    elseif i >= 20000 && i <= 30000
        F = 100;
    elseif i >= 40000 && i <= 50000
        F = 100;
    elseif i >= 60000 && i <= 70000
        F = 100;
    else
        F = 0;
    end
    t = Kod*i; %Scale t here!

    %ODES: Michaelis-Menten
    k1C = ((1+F)/(1+(o(i)/Ki1)))-kcd*c(i);
    k1A = ((c(i))/(1+(o(i)*r(i)/Ki2)))-kad*a(i);
    k1R = (((o(i)*r(i))^2)/(K+(o(i)*r(i))^2)) + Kor - krd*r(i);
    k1O = a(i) - o(i);
    
    k2C = (((1+F)/(1+(o(i)+0.5*k1O*dt)/Ki1)))-kcd*(c(i)+0.5*k1C*dt);
    k2A = (((c(i)+0.5*k1C*dt))/(1+(o(i)+0.5*k1O*dt)*r(i)/Ki2))-kad*(a(i)+0.5*k1A*dt);
    k2R = (((o(i)+0.5*k1O*dt)*(r(i)+0.5*k1R*dt))^2)/(K+(o(i)+0.5*k1O*dt)*(r(i)+0.5*k1R*dt)^2) + Kor - krd*(r(i)+0.5*k1R*dt);
    k2O = (a(i)+0.5*k1A*dt) - (o(i)+0.5*k1O*dt);
    
    k3C = (((1+F)/(1+(o(i)+0.5*k2O*dt)/Ki1)))-kcd*(c(i)+0.5*k2C*dt);
    k3A = (((c(i)+0.5*k2C*dt))/(1+(o(i)+0.5*k2O*dt)*r(i)/Ki2))-kad*(a(i)+0.5*k2A*dt);
    k3R = (((o(i)+0.5*k2O*dt)*(r(i)+0.5*k2R*dt))^2)/(K+(o(i)+0.5*k2O*dt)*(r(i)+0.5*k2R*dt)^2) + Kor - krd*(r(i)+0.5*k2R*dt);
    k3O = (a(i)+0.5*k2A*dt) - (o(i)+0.5*k2O*dt);
    
    k4C = (((1+F)/(1+(o(i)+k3O*dt)/Ki1)))-kcd*(c(i)+k3C*dt);
    k4A = (((c(i)+k3C*dt))/(1+(o(i)+k3O*dt)*r(i)/Ki2))-kad*(a(i)+k3A*dt);
    k4R = (((o(i)+k3O*dt)*(r(i)+k3R*dt))^2)/(K+(o(i)+k3O*dt)*(r(i)+k3R*dt)^2) + Kor - krd*(r(i)+k3R*dt);
    k4O = (a(i)+k3A*dt) - (o(i)+k3O*dt);
    
    c(i+1) = c(i) + (1/6)*(k1C+2*k2C+2*k3C+k4C)*dt;
    a(i+1) = a(i) + (1/6)*(k1A+2*k2A+2*k3A+k4A)*dt;
    o(i+1) = o(i) + (1/6)*(k1O+2*k2O+2*k3O+k4O)*dt;
    r(i+1) = r(i) + (1/6)*(k1R+2*k2R+2*k3R+k4R)*dt;
    T(i+1) = T(i) + dt;
end
% Make 5 different plots charts for c, a, o, r in one figure
figure(1)
subplot(5,1,1)
plot(T,o)
title('Cortisol')
xlabel('Time (Hours)')
ylabel('concentration')

subplot(5,1,2)
plot(T,a)
title('Adrenocorticotropic hormone')
xlabel('Time (Hours)')
ylabel('concentration')

subplot(5,1,3)
plot(T,c)
title('Corticotropin releasing hormone')
xlabel('Time (Hours)')
ylabel('concentration')

subplot(5,1,4)
plot(T,r)
title('Glucocorticoid Receptor')
xlabel('Time (Hours)')
ylabel('concentration')

