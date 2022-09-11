function EEsteadyStateforSVEIRS
%steadyState symbolic calculation for the steady state of an SVEIIRS model


    disp('SVEIRS model')    
    syms S V E I0 I1 I2 I3 R N nu b0 b1 b2 b3 N0 taue tau0 tau1 tau2 tau3 R0 c c0 c1 c2 c3 Nc etar etav


    N = S+ V+E+ I0 + I1 + I2 + I3 +R % = N0
    Nc = S*c + V*c + +E*c +I0*c0 + I1*c1+ I2*c2+ I3*c3+ R*c
    alphas0 = b0*c0*c*S/Nc
    alphas1 = b1*c1*c*S/Nc
    alphas2 = b2*c2*c*S/Nc
    alphas3 = b3*c3*c*S/Nc
    
    % dSdt = etar*R + etav*V -alpha0*I_0 -alpha1*I_1 -alpha2*I_2 -alpha3*I_3-nu*S;
    % dVdt = nu*S - etav*V;
    % dEdt = alpha0*I_0 +alpha1*I_1 +alpha2*I_2 +alpha3*I_3 - gammae*E;
    % dI_0dt = gammae*E-gamma0*I_0;
    % dI_1dt = gamma0*I_0-gamma1*I_1;
    % dI_2dt = gamma1*I_1-gamma2*I_2;
    % dI_3dt = gamma2*I_2-gamma3*I_3;
    % dRdt = gamma3*I_3 - etar*R;
    % dTdt = alpha0*I_0 +alpha1*I_1 +alpha2*I_2 +alpha3*I_3;

A = [ 0,   0, 1, 0, 0, 0, 0,  -taue*etar;
      0,   0, 0, 1, 0, 0, 0,  -tau0*etar;
      0,   0, 0, 0, 1, 0, 0,  -tau1*etar;
      0,   0, 0, 0, 0, 1, 0,  -tau2*etar;
      0,   0, 0, 0, 0, 0, 1,  -tau3*etar;
      1 - 1/R0, -(1/R0), 0, 0, 0, 0, 0, -((etar*taue)/R0 + (etar*(c0*tau0 + c1*tau1 +c2*tau2 +c3*tau3))/(c*R0) + 1/R0);
      1, -etav/nu, 0, 0, 0, 0, 0, 0;
      1, 1, 1, 1, 1, 1, 1, 1] 
 

B = [0; 
    0; 
    0; 
    0;  
    0;
    0;
    0;
    N0]

X = linsolve(A,B)
S = X(1)
V = X(2)
E = X(3)
I0 = X(4)
I1 = X(5)
I2 = X(6)
I3 = X(7)
R = X(8)
N = S+ V+ E+ I0 + I1 + I2 + I3 +R


end