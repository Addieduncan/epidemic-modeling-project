function next_generation_matrix_R0forSVEIRS
%next_gen_matrix_R0 symbolic calculation for the Next Generation Matrix of
%an SVEIRS model

% In the next generation matrix approach, we define the infection vector
% X = (I1, I2, ...)'
% where I* are the infected components. and express the ODEs as
% dX/dt = F(X) - V(X)
% where
%
% F(X)= vector function of all the new infections entering the system
% V(X)= vector funtion of all the progressions between the components
%
% We then define
%
% J_f = Jacobian of the new infection function
% J_v = Jacobian of the progression function
%
% After linearizing the equations about the disease free equilibrium (DFE)
% where all the infection states are zero, X = (0, 0, ...), the linearized
% system of equations can be written as
% dX/dt = J_f*X  - J_v*X  for X ~ 0
%       = (J_f*J_v^{-1}  - 1)*J_v*X
%
% The next generation matrix is defined as NGM = J_f*J_v^{-1} and this
% system of ODEs is stable if all the eigenvalues of NGM < 1.  The largest
% eigenvalue of the NGM is the basic reproductive number Ro.
% independent variables
% alphaI = force of infection from infectious compartment I
% tauI = mean time spent in compartment I
% tauM = migration time
% PsI = probability of progressing to infected state I
%
% dependent variables
%  gammaI = 1/tau* progression rate for leaving compartment * if tau_m=0
% Pij = probability that the person in compartment i progresses to compartment j
% lambda = force of infection on the susceptible population

%  September 21, 2020, Mac



    disp('SVEIRS model')
    
    syms E I0 I1 I2 I3 alpha0 alpha1 alpha2 alpha3 taue tau0 tau1 tau2 tau3 gamma0 gamma1 gamma2 gamma3 gammae 
    
    mu=0; % no migration
    X=[E;I0;I1;I2;I3]
    
    F=[alpha0*I0 + alpha1*I1 + alpha2*I2 + alpha3*I3;
        0;
        0;
        0;
        0]
    
    V = [ E*gammae
        -gammae*E + I0*gamma0;
        -I0*gamma0 + I1*gamma1;
        -I1*gamma1 + I2*gamma2;
        -I2*gamma2 + I3*gamma3]
    

J_f = jacobian(F,X)
J_v = jacobian(V,X)

J_v_inv = inv(J_v)

NGM = J_f*J_v_inv
N_ev = eig(NGM)

end