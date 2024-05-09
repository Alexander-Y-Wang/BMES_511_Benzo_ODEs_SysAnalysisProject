%% RK4 MODE
clear all 
%% Initial 
%% Initial 
C = [120]; %corticotropin releasing hormone
A = [1.2]; %adrenocorticotropic hormone
O = [0.0065]; %cortisol
R = [.08]; %glucocorticoid receptor on adrenal
OR = [.1]; %homodimer of GR bound to cortisol (O) *Stays in quasi-equilibrium*

t = [0]; %initial time
dt = 0.0001; %time step
t_end = 50; %time end

%% Constants
%CRK (C): Corticotropin releasing hormone
Kc = 1; %Production of C
Kcd = 1; %Constant degradation of C
Kor = 0.05; %pituitary GR production

% a (ACTH): adrenocorticotropic hormone
Ka = 1; %Production of a
Kad = 10; %Constant degradation of a

% o (O): Cortisol
Ko = 1; %Production of o
Kod = 1; %Constant degradation of o

% r (R): glucocorticoid receptor (IN PITIUTARY) 
Kr = 1; %Production of r
Krd = 0.9;  %Constant degradation of r

% i : inhibition constants
Ki1 = 0.1; %inhibition 1
Ki2 = 0.1; %inhibition 2

K = 0.001; %equilibrium binding affinity


for i = 1:(t_end/dt)
    if i <= 100000
        F = 7600; %This stress value was pulled from Srisam (PTSD paper)
    elseif i >= 300000 && i <= 400000
        F = 7600;
    else
        F = 0;
    end
    k1c = (Kc+F)/(1+(O(i)/Ki1))-Kcd*C(i);
    k2c = (Kc+F)/(1+(O(i)+0.5*k1c*dt)/Ki1)-Kcd*(C(i)+0.5*k1c*dt);
    k3c = (Kc+F)/(1+(O(i)+0.5*k2c*dt)/Ki1)-Kcd*(C(i)+0.5*k2c*dt);
    k4c = (Kc+F)/(1+(O(i)+k3c*dt)/Ki1)-Kcd*(C(i)+k3c*dt);
    C(i+1) = C(i) + (dt/6)*(k1c + 2*k2c + 2*k3c +k4c);
    
    k1a = (Ka*C(i))/(1+OR(i)/Ki2)-Kad*A(i);
    k2a = (Ka*(C(i)+0.5*k1c*dt))/(1+(OR(i)+0.5*k1a*dt)/Ki2)-Kad*(A(i)+0.5*k1a*dt);
    k3a = (Ka*(C(i)+0.5*k2c*dt))/(1+(OR(i)+0.5*k2a*dt)/Ki2)-Kad*(A(i)+0.5*k2a*dt);
    k4a = (Ka*(C(i)+k3c*dt))/(1+(OR(i)+k3a*dt)/Ki2)-Kad*(A(i)+k3a*dt);
    A(i+1) = A(i) + (dt/6)*(k1a + 2*k2a + 2*k3a +k4a);
    
    k1o = Ko*A(i) - Kod*O(i);
    k2o = Ko*(A(i)+0.5*k1a*dt) - Kod*(O(i)+0.5*k1o*dt);
    k3o = Ko*(A(i)+0.5*k2a*dt) - Kod*(O(i)+0.5*k2o*dt);
    k4o = Ko*(A(i)+k3a*dt) - Kod*(O(i)+k3o*dt);
    O(i+1) = O(i) + (dt/6)*(k1o + 2*k2o + 2*k3o +k4o);
    
    k1r = (Kr*(OR(i))^2)/(K+(OR(i))^2) + Kor - Krd*R(i);
    k2r = (Kr*(OR(i)+0.5*k1r*dt)^2)/(K+(OR(i)+0.5*k1r*dt)^2) + Kor - Krd*(R(i)+0.5*k1r*dt);
    k3r = (Kr*(OR(i)+0.5*k2r*dt)^2)/(K+(OR(i)+0.5*k2r*dt)^2) + Kor - Krd*(R(i)+0.5*k2r*dt);
    k4r = (Kr*(OR(i)+k3r*dt)^2)/(K+(OR(i)+k3r*dt)^2) + Kor - Krd*(R(i)+k3r*dt);
    R(i+1) = R(i) + (dt/6)*(k1r + 2*k2r + 2*k3r +k4r);
    OR(i+1) = OR(i); %*Stays in quasi-equilibrium, thus no change is appropriate*
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