%% RK4 MODE
%% Initial 

C = [120]; %corticotropin releasing hormone
A = [1.2]; %adrenocorticotropic hormone
O = [0.0065]; %cortisol
R = [1.065]; %glucocorticoid receptors on adrenal
OR = [.1]; %homodimer of GR bound to cortisol (O) *Stays in quasi-equilibrium*

t = [0]; %initial time
dt = 0.0001; %time step
t_end = 50; %time end

%% Constants
%CRK (C): Corticotropin releasing hormone
Kc = 1; %Production of Corticotropin releasing hormone
Kcd = 10; %Constant degradation of Corticotropin releasing hormone
Kor = 0.05; %pituitary GR production

% a (ACTH): adrenocorticotropic hormone
Ka = 1; %Production of adrenocorticotropic hormone
Kad = 10; %Constant degradation of adrenocorticotropic hormone

% o (O): Cortisol
Ko = 1; %Production of Cortisol
Kod = 1; %Constant degradation of Cortisol

% r (R): glucocorticoid receptor (IN PITIUTARY) 
Kr = 1; %Production of glucocorticoid receptor
Krd = .9;  %Constant degradation of glucocorticoid receptor !!!!!!!!!!!!

% i : inhibition constants
Ki1 = 0.1; %inhibition 1
Ki2 = 0.1; %inhibition 2

K = 0.001; %equilibrium binding affinity

for i = 1:(t_end/dt)
    if i <= 100000
        F = 7600; %This stress value was pulled from Srisam (PTSD paper)
    elseif i >= 300000 && i <= 400000 %Simulate some applied stresses
        F = 7600;
    else
        F = 0;
    end

    %ODES: Michaelis-Menten
    k1C = (Kc+F)/(1+(O(i)/Ki1))-Kcd*C(i);
    k1A = (Ka*C(i))/(1+(OR(i)/Ki2))-Kad*A(i);
    k1R = (Kr*(OR(i)^2))/(K+(OR(i)^2)) + Kor - Krd*R(i);
    k1O = Ko*A(i) - Kod*O(i);
    
    k2C = ((Kc+F)/(1+(O(i)+0.5*k1O*dt)/Ki1))-Kcd*(C(i)+0.5*k1C*dt);
    k2A = ((Ka*(C(i)+0.5*k1C*dt))/(1+(OR(i)+0.5*k1R*dt)/Ki2))-Kad*(A(i)+0.5*k1A*dt);
    k2R = (Kr*((OR(i)+0.5*k1R*dt)^2))/(K+((OR(i)+0.5*k1R*dt)^2)) + Kor - Krd*(R(i)+0.5*k1R*dt);
    k2O = Ko*(A(i)+0.5*k1A*dt) - Kod*(O(i)+0.5*k1O*dt);
    
    k3C = ((Kc+F)/(1+(O(i)+0.5*k2O*dt)/Ki1))-Kcd*(C(i)+0.5*k2C*dt);
    k3A = ((Ka*(C(i)+0.5*k2C*dt))/(1+(OR(i)+0.5*k2R*dt)/Ki2))-Kad*(A(i)+0.5*k2A*dt);
    k3R = (Kr*((OR(i)+0.5*k2R*dt)^2))/(K+((OR(i)+0.5*k2R*dt)^2)) + Kor - Krd*(R(i)+0.5*k2R*dt);
    k3O = Ko*(A(i)+0.5*k2A*dt) - Kod*(O(i)+0.5*k2O*dt);
    
    k4C = ((Kc+F)/(1+(O(i)+k3O*dt)/Ki1))-Kcd*(C(i)+k3C*dt);
    k4A = ((Ka*(C(i)+k3C*dt))/(1+(OR(i)+k3R*dt)/Ki2))-Kad*(A(i)+k3A*dt);
    k4R = (Kr*((OR(i)+k3R*dt)^2))/(K+((OR(i)+k3R*dt)^2)) + Kor - Krd*(R(i)+k3R*dt);
    k4O = Ko*(A(i)+k3A*dt) - Kod*(O(i)+k3O*dt);
    
    C(i+1) = C(i) + (1/6)*(k1C+2*k2C+2*k3C+k4C)*dt;
    A(i+1) = A(i) + (1/6)*(k1A+2*k2A+2*k3A+k4A)*dt;
    O(i+1) = O(i) + (1/6)*(k1O+2*k2O+2*k3O+k4O)*dt;
    R(i+1) = R(i) + (1/6)*(k1R+2*k2R+2*k3R+k4R)*dt;
    OR(i+1) = OR(i); %*Stays in quasi-equilibrium*
    t(i+1) = t(i) + dt;
end


% Make 5 different plots charts for C, A, O, R, and OR in one figure
figure
subplot(5,1,1)
plot(t,C)
title('Corticotropin releasing hormone')
xlabel('Time (Hours)')
ylabel('concentration')

subplot(5,1,2)
plot(t,A)
title('Adrenocorticotropic hormone')
xlabel('Time (Hours)')
ylabel('concentration')

subplot(5,1,3)
plot(t,O)
title('Cortisol')
xlabel('Time (Hours)')
ylabel('concentration')

subplot(5,1,4)
plot(t,R)
title('Glucocorticoid receptor on pituitary')
xlabel('Time (Hours)')
ylabel('GR')

subplot(5,1,5)
plot(t,OR)
title('homodimer of GR bound to cortisol')
xlabel('Time (Hours)')
ylabel('OR')