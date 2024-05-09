clear all
%Model code, trying to figure out how to Euler it. Variables are subject
%to change
%% EULER MODE
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
Kcd = 10; %Constant degradation of C
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


    dCdt = F*(K)


    C(i+1) = C(i) + dt*dCdt;
    A(i+1) = A(i) + dt*dAdt;
    O(i+1) = O(i) + dt*dOdt;
    R(i+1) = R(i) + dt*dRdt;
    OR(i+1) = OR(i); %*Stays in quasi-equilibrium*
    t(i+1) = t(i) + dt;
    
end

% Plot results
% Make 5 different plots charts for C, A, O, R, and OR in one figure
figure
subplot(5,1,1)
plot(t,C)
title('Corticotropin releasing hormone')
xlabel('Time (Hours)')
ylabel('concentration (pmol/L)')

subplot(5,1,2)
plot(t,A)
title('Adrenocorticotropic hormone')
xlabel('Time (Hours)')
ylabel('concentration (pmol/L)')

subplot(5,1,3)
plot(t,O)
title('Cortisol')
xlabel('Time (Hours)')
ylabel('concentration (pmol/L)')

subplot(5,1,4)
plot(t,R)
title('Glucocorticoid receptor on pituitary')
xlabel('Time (Hours)')
ylabel('GR')

