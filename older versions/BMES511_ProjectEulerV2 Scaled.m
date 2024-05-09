clear all
%% EULER MODE OF SCALED VERSION!!!! <-----------------
%% Initial 

C(1) = 120; %corticotropin releasing hormone
A(1) = 1.2*0.001; %adrenocorticotropic hormone
O(1) = 6*0.001; %cortisol
R(1) = .1; %glucocorticoid receptor on adrenal
OR(1) = .1; %homodimer of GR bound to cortisol (O) *Stays in quasi-equilibrium*

t(1) = 0; %initial time
dt = 0.1; %time step
t_end = 120; %time end

%% Constants
%CRK (C): Corticotropin releasing hormone
Kc = 1; %Production of C
Kcd = 1; %Constant degradation of C
Kcr = 0.05; %pituitary GR production

% a (ACTH): adrenocorticotropic hormone
Ka = 10; %Production of a
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



%Scaled variables
C(1) = Kod*C(1)/Kc; %corticotropin releasing hormone
A(1) = (Kod^2)*A(1)/(Kc*Ka); %adrenocorticotropic hormone
O(1) = ((Kod^3)*O(1))/(Kc*Ka*Ko); %cortisol
R(1) = (Kod*R(1))/Kr; %glucocorticoid receptor on adrenal
OR(1) = .1; %homodimer of GR bound to cortisol (O) *Stays in quasi-equilibrium*


for i = 1:(t_end/dt)
    F = 0.02 + 0.01*sin(2*pi*t(i)/24); % external stress stimuli (sinusoidal function)
    
    %Scaled t:
    t = Kod*t(i)

    %ODES: Michaelis-Menten SCALED 
    dCdt = (1+F)/(1+(O(i)/Ki1))-Kcd*C(i); %hypothalumus
    dAdt = C(i)/(1+OR(i)/Ki2)-Kad*A(i); %pituitary
    dRdt = OR(i)^2)./(K+(OR(i)).^2) + Kcr - Krd*R(i); %pituitary
    dOdt = A(i) - Kod*O(i); %adrenal
    
    C(i+1) = C(i) + dt*dCdt;
    A(i+1) = A(i) + dt*dAdt;
    O(i+1) = O(i) + dt*dOdt;
    R(i+1) = R(i) + dt*dRdt;
    OR(i+1) = OR(i); %*Stays in quasi-equilibrium*
    t(i+1) = t(i) + dt;
end

plot(t,C,t,A,t,O,t,R,t,OR)
legend('C','A','O','R','OR')