function SVEIRSmodel
%SIRepidemicModel - solves and SIR model 

%   Detailed explanation goes here


% set default values for graphics and printing
clear global ; clf; format shortE; close all;  % close previous sessions
set(0,'DefaultAxesFontSize',18);set(gca,'FontSize',18);close(gcf); % increase font size
%linespec = {'-b','-.g','--r',':c','-m*','-.bo','--g+',':rs','-cd'};
linespec = {'-b','-.g','--r',':c',':m',':b',':g','-.m','-c'};

% define the independent model parameters
c = 11 ; % contacts per day for S, V, E, R
c0 = 11 ; % contacts per day for I0
c1 = 3 ; % contacts per day for I1
c2 = 7 ; % contacts per day for I2
c3 = 11 ; % contacts per day for I3

rho = .7; % mask coefficient


b0 = 0.02; % base infectivity of I0
b1= 0.09; % base infectivity of I1
b2= 0.05; % base infectivity of I2
b3= 0.02; % base infectivity of I3
beta0 = b0*(1-rho); % infectivity of I0 in a pandemic
beta1= b1*(1-rho); % infectivity of I1 in a pandemic
beta2= b2*(1-rho); % infectivity of I2 in a pandemic
beta3= b3*(1-rho); % infectivity of I3 in a pandemic
%beta0 = b0; % infectivity of I0
%beta1= b1; % infectivity of I1
%beta2= b2*(1-rho); % infectivity of I2
%beta3= b3*(1-rho); % infectivity of I3

taue = 2; % average time a person is latent
tau0 = 1; % average time a person is in I0
tau1 = 2; % average time a person is in I1
tau2 = 2; % average time a person is in I2
tau3 = 3; % average time a person is in I3

nu = 0; % vaccination rate

zetar = 240; % average time a recovered person is immune
zetav = 240; %average time a vaccinated person is immune

N0 = 1; % total population (normalized)
E0 = 0.0020; % initial population
I_00 = 0.0020; % initial population
I_10 = 0.0020; % initial population
I_20 = 0.0020; % initial population
I_30 = 0.0020; % initial population
T0 = E0+ I_00 + I_10 + I_20 + I_30; % initial incidence
V0=0; % initial population
R0=0; % initial population
S0=N0- T0 - V0 - R0;% initial susceptible population

varlabels = {'S','V','E','I_0','I_1','I_2','I_3','R','T'};
nvar=length(varlabels);

% define the independent simulation parameters
tfinal= 365*3; % final time in days
nout = 365*3+1; % number of times the output is saved

% define the dependent model parameters
gammae = 1/taue; % recovery rate for E
gamma0 = 1/tau0; % recovery rate for I0
gamma1 = 1/tau1; % recovery rate for I1
gamma2 = 1/tau2; % recovery rate for I2
gamma3 = 1/tau3; % recovery rate for I3

etar = 1/zetar; % loss of immunity rate for R
etav = 1/zetav; % loss of immunity rate for R

Ps = 1; %probability that a contact is susceptible at being of epidemic

alpha0=c0*beta0*Ps; % force from infection for I0
alpha1=c1*beta1*Ps; % force from infection for I1
alpha2=c2*beta2*Ps; % force from infection for I2
alpha3=c3*beta3*Ps; % force from infection for I3



% define the dependent simulation parameters
tspan = linspace(0,tfinal,nout);
y0=copyvar(S0,V0,E0,I_00,I_10,I_20,I_30,R0,T0); % initial conditions for the ODE solver

% analyze problem and initial conditions
R0 = alpha0/gamma0 +alpha1/gamma1 + alpha2/gamma2 + alpha3/gamma3; % formula for the basic reproductive number
disp(['R0 = ',num2str(R0)])

%
ode = @(t,y) functionODE(t,y,c,c0,c1,c2,c3,beta0,beta1,beta2,beta3,gamma0,gamma1,gamma2,gamma3,gammae,etar,etav,nu,rho,N0);
[tout,yout] = ode45(ode, tspan, y0); % ode45 ode113
disp([' tout ',varlabels])
disp([tout,yout])

% plot the solution
for ip=1:nvar
plot(tout,yout(:,ip),linespec{ip}); hold on;
end
xlabel('time')
ylabel('populations')
title(['R_0 = ',num2str(R0)])
legend(varlabels)

% compute the time derivatives
dydttout=NaN(size(yout));
for it=1:nout
    dydttout(it,:)=functionODE(tout(it),yout(it,:),c,c0,c1,c2,c3,beta0,beta1,beta2,beta3,gamma0,gamma1,gamma2,gamma3,gammae,etar,etav,nu,rho,N0);
end

% plot the derivatives of the solution
figure
for ip=1:nvar
plot(tout,dydttout(:,ip),linespec{ip}); hold on;
end
xlabel('time')
ylabel('derivatives')
title(['R_0 = ',num2str(R0)])
legend(varlabels)

end

function dydt=functionODE(t,y,c,c0,c1,c2,c3,beta0,beta1,beta2,beta3,gamma0,gamma1,gamma2,gamma3,gammae,etar,etav,nu,rho,N0)
%functionODE define the time derivative

% extract the physically meaningful variables
[S, V,E, I_0, I_1, I_2, I_3, R , T]=copyvar(y);

% probability that a contact is susceptible
Ps = c*S/(c*S + c*V + c*E +c0*I_0 +c1*I_1 +c2*I_2 +c3*I_3 + c*R); % nonlinear model

alpha0=c0*beta0*Ps; % force from infection for I0
alpha1=c1*beta1*Ps; % force from infection for I1
alpha2=c2*beta2*Ps; % force from infection for I2
alpha3=c3*beta3*Ps; % force from infection for I3


dSdt = etar*R + etav*V -alpha0*I_0 -alpha1*I_1 -alpha2*I_2 -alpha3*I_3-nu*S;
dVdt = nu*S - etav*V;
dEdt = alpha0*I_0 +alpha1*I_1 +alpha2*I_2 +alpha3*I_3 - gammae*E;
dI_0dt = gammae*E-gamma0*I_0;
dI_1dt = gamma0*I_0-gamma1*I_1;
dI_2dt = gamma1*I_1-gamma2*I_2;
dI_3dt = gamma2*I_2-gamma3*I_3;
dRdt = gamma3*I_3 - etar*R;
dTdt = alpha0*I_0 +alpha1*I_1 +alpha2*I_2 +alpha3*I_3;

% define the time derivatives for the ODE solver
dydt=copyvar(dSdt,dVdt,dEdt,dI_0dt,dI_1dt,dI_2dt,dI_3dt,dRdt,dTdt);
end
function varargout=copyvar(varargin)
%copyvar copies varargin to varargout
% the input or output can be an array or list of variables

% 3 cases to consider
if nargout == nargin % same number output as input
    varargout=varargin;
    
elseif nargout < nargin % output is an nvarin x 1 array
    y=NaN(nargin,1);
    for i=1:nargin
        y(i)=varargin{i};
    end
    varargout{1}=y;
    
else                     % input is an nvarin x 1 array
    yin=varargin{1};
    nyin=length(yin);
    varargout=cell(nargout,1);
    for i=1:nyin
        varargout{i}=yin(i);
    end
    
end

end