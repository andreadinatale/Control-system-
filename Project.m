%% Project Aerospace Control System A.Y 2020/2021

clear all;
close all;
clc;

%% DATA

s=tf('s');
g=9.81;

% Stability derivatives

Yv=ureal('Yv',-0.264,'Perc',3*4.837);
Yp=0;
Lv=ureal('Lv',-7.349,'Perc',3*4.927);
Lp=0;

% Control derivatives

Yd=ureal('Yd',9.568,'Perc',3*4.647);
Ld=ureal('Ld',1079.339,'Perc',3*2.762);

%  N.B====> use 3*sigma in the Gaussian distribution because our
%  requirements must be satisfied in +-3*sigma

%% DEFINE THE LATERAL DYNAMICS STATE SPACE SYSTEM

A=[Yv Yp g
    Lv Lp 0
    0 1 0];
B=[Yd Ld 0]';
C=[0 1 0
    0 0 1];
D=[0 0]';
states={'v','p','phi'};
input={'Delta lat'};
outputs={'p','phi'};

G=ss(A,B,C,D,'statename',states,'inputname',input,'outputname',outputs);  % uncertain plant in SS form
Gn=G.NominalValue;   % nominal plant in SS form

figure('Name','Lateral dynamics','NumberTitle','off'), bode(G), grid on, hold on;
bode(Gn,'r'), legend('Uncertain plant','Nominal plant');
title('Lateral dynamics');
figure('Name','Step response of uncertain and nominal plan','NumberTitle','off'), step(G), grid on, hold on 
step(Gn), legend('Uncertain plant','Nominal plant')   
title('Step response of uncertain and nominal plant');

disp('PROPERTIES OF LATERAL DYNAMICS');
damp(G)   % is unstable!!!!!

figure('Name','Poles and zeros of uncertain and nominal plant','NumberTitle','off'), pzmap(G), grid on, hold on;
pzmap(Gn), legend('Uncertain plant','Nominal plant','Location','East');
title('Poles and zeros of uncertain and nominal plant');

%% DESIGN APPROACHES: SYSTUNE 

% Define tunable controllers

R_phi_tun=tunablePID('R_phi_P','P');
R_phi_tun.InputName='ephi';            % tunable R_phi
R_phi_tun.OutputName='p0';

R_p_tun=tunablePID2('R_p_PID','PID');
R_p_tun.InputName={'p0','p'};                  % tunable R_p
R_p_tun.OutputName='Delta lat';

% Build the closed loop system from phi0 to [p phi]'

sum_ephi=sumblk('ephi=phi0-phi'); % create the summation junction of error on phi

CL0=minreal(connect(Gn,R_p_tun,R_phi_tun,sum_ephi,'phi0',{'phi','Delta lat','ephi'}));  % generalized continuous-time state-space model from phi0 to [phi Delta lat ephi]'

% Specify the requirements: 
% 1)second order step responce with omegan>=10 rad/s and xsi>=0.9
% 2)Control effort limitation: |Delta lat|<5 deg for a doublet step of phi0

omegan=10;
xsi=0.9;

ref_sys=1/(1+2*(xsi/omegan)*s+(s/omegan)^2);  % Response requirements: second order transfert function 
Wp_inv=1-ref_sys;
Wp=1/Wp_inv;
Wq_inv=tf([0.25 20*0.25],[1 10]); % Control effort limitaition: |Delta lat|<5 deg for a doublet step of phi0
figure('Name','Wq_inv','NumberTitle','off'), bodemag(Wq_inv), grid on, legend('Wq_inv');

requirements=[TuningGoal.StepTracking('phi0','phi',ref_sys) TuningGoal.WeightedGain('phi0','Delta lat',Wq_inv,[])];  % construct the requirements 

[CL,fSoft,~,info_CL]=systune(CL0,requirements);   % tune the closed loop system 

figure('Name','Step response','NumberTitle','off'), viewGoal(requirements,CL);

%% DEFINED THE TUNED CONTROLLERS

disp('VALUES OF TUNED PID CONTROLLERS');
Rphi=pid(CL.Blocks.R_phi_P.Kp.Value,CL.Blocks.R_phi_P.Ki.Value,CL.Blocks.R_phi_P.Kd.Value,CL.Blocks.R_phi_P.Tf.Value)
Rphi.InputName='ephi';
Rphi.OutputName='p0';
Rp=pid2(CL.Blocks.R_p_PID.Kp.Value,CL.Blocks.R_p_PID.Ki.Value,CL.Blocks.R_p_PID.Kd.Value,CL.Blocks.R_p_PID.Tf.Value,CL.Blocks.R_p_PID.b.Value,CL.Blocks.R_p_PID.c.Value)
Rp.InputName={'p0','p'};
Rp.OutputName='Delta lat';

%% DEFINE THE UNCERTAIN CLOSED LOOP TF

CL_unce=minreal(connect(G,Rp,Rphi,sum_ephi,'phi0',{'phi','Delta lat','ephi'}));

%% 

figure('Name','Nominal design: Sensitivity weight','NumberTitle','off'), bodemag(CL_unce(3)), grid on, hold on;
bodemag(Wp_inv,'--r'), legend('Sensitivity','Wp inv');
title('Nominal design: Sensitivity weight');
legend('Location','Southeast');

figure('Name','Nominal design: Control sensitivity','NumberTitle','off'), bodemag(CL_unce(2)), grid on, hold on;
bodemag(Wq_inv); 
title('Nominal design: Control sensitivity');
legend('Control effort','Wq inv');

figure('Name','Step response','NumberTitle','off'), step(CL_unce), grid on;

%% CHECK ROBUST STABILITY (RS)

% Graphically, I need to verify that Fn<1/W, with Wrs>=|G-Gn|/Gn

% Compute relative error

Garray=usample(G(2),70);
[P,info_RS]=ucover(Garray,Gn(2),4);  % create a 4° order W
W=info_RS.W1;

% Verify RS graphically

figure('Name','Relative error and weight for Robust stability','NumberTitle','off'), bodemag((Garray-Gn(2))/Gn(2),W,'g'), grid
legend('Relative error','W');
title('Relative error and weight for Robust stability');

figure('Name','Robust stability verification','NumberTitle','off'), bodemag(CL(1),'r',1/W,'g',[10^(-5):10^3]), grid on;
legend('Complementary sensitivity','1/W');
title('Robust stability verification');

% RS of nominal plan is guaranteed!!!

%% CHECK THE PERFORMACE RELATED TO DELTA LAT BOUNDED FROM A DOUBLET INPUT OF PHI0

% Define the input of phi0 both in time and Laplace domain

phi0_t=@(t) (0.*t).*(t<=1)+(10+0.*t).*(t>=1).*(t<=3)+(-10+0.*t).*(t>=3).*(t<=5)+(0.*t).*(t>=5);  % NB [deg]
t=[0:0.001:6];
figure('Name','Phi0 input','NumberTitle','off'), plot(t,phi0_t(t)), grid on;
xlabel('Time (seconds)');
ylabel('Phi0 (rad)');
title('Phi0 input');

figure('Name','Doublet Response','NumberTitle','off'), lsim(CL_unce,phi0_t(t),t), grid on

%% REPRESENT THE UNCERTAIN CLOSED LOOP SYSTEM BOTH FOR SYSTUNE

% Define the uncertain closed loop system from phi0 to [p phi]'

figure('Name','Complementary sensitivity','NumberTitle','off'), bodemag(CL_unce(1)), grid on;
title('Complementary sensitivity');

figure('Name','Step response of uncertain model','NumberTitle','off'), step(CL_unce(1)), grid on, hold on;
step(ref_sys), legend('Actual','Desired');
title('Step response of uncertain model');

%% ROBUST STABILITY AND PERFORMANCE VERIFICATION     

optrp=robOptions('Display','on');
[perfmarg,wcu,info]=robgain(CL_unce(1),1.01,optrp)  % calculate the performance margin in terms of lower and upper bound
                                                   % and wcu is the worst
                                                   % combination of the
                                                   % uncertain parameters
%% M-DELTA DECOMPOSITION

% omega=logspace(-4,5,1500);
%
[M,Delta,Blocks]=lftdata(CL_unce);
szDelta=size(Delta);
M_red=M(1:szDelta(2),1:szDelta(1));
% M_red_w=frd(M_red,omega);

LinMagopt = bodeoptions;
LinMagopt.PhaseVisible = 'off'; LinMagopt.XLim = [1e-1 1e2]; LinMagopt.MagUnits = 'abs';
mu_bounds=mussv(M_red,Blocks,'s');
figure('Name','Structured singolar value','NumberTitle','off'), bodeplot(mu_bounds(1,1),mu_bounds(1,2),LinMagopt), grid on;
xlabel('Frequency');
ylabel('Mu upper/lower bounds');
title('Structured singolar value');
legend('Mu upper bound','Mu lower bound');

figure('Name','Singolar values of M','NumberTitle','off'), sigma(mu_bounds), grid on
title('Singolar values of M');

%% MONTECARLO SIMULATION

N=1000;

F_st_mc=zeros(N,1);
F_rt_mc=zeros(N,1);
F_os_mc=zeros(N,1);
delta_lat_mc=zeros(N,1);
Gm=zeros(N,1);
Pm=zeros(N,1);
Sm=zeros(N,1);

for n=1:N
    % Complementary sensitivity
    
    F_0=usample(CL_unce(1),1);
    
    [y_F,t_F]=step(F_0,linspace(0,10,1000));
    
    step_info_F=stepinfo(y_F,t_F,1);
    
    F_st_mc(n)=step_info_F.SettlingTime;
    F_rt_mc(n)=step_info_F.RiseTime;
    F_os_mc(n)=step_info_F.Overshoot;
    
    % Control effort
    
    Q_0=usample(CL_unce(2),1);
    
    [y_Q,t_Q]=lsim(Q_0,phi0_t(linspace(0,10,1000)),linspace(0,10,1000));
    
    delta_lat_mc(n)=max(abs(y_Q));
    
    % Loop transfer function L=F/(1-F)
    
    L=CL_unce(1)/(1-CL_unce(1));
    S=1/(1+L);
    S_0=usample(S,1);
    [Gm(n),Pm(n)]=margin(F_0/S_0);
    Sm(n)=1/getPeakGain(S_0);
    
end

figure('Name','Montecarlo: Phase margin','NumberTitle','off'), histogram(Pm,100), grid on, title('Phase margin');
xlabel('Phase margin');
figure('Name','Montecarlo: Stability margin','NumberTitle','off'), histogram(Sm,100), grid on, title('Stability margin');
xlabel('Stability margin');
figure('Name','Montecarlo: Settling time','NumberTitle','off'), histogram(F_st_mc,100), grid on, title('Settling time');
xlabel('Settling time');
figure('Name','Montecarlo: % Overshoot','NumberTitle','off'), histogram(F_os_mc,100), grid on, title('% Overshoot');
xlabel('% Overshoot');
figure('Name','Montecarlo: Rise Time','NumberTitle','off'), histogram(F_rt_mc,100), grid on, title('Rise Time');
xlabel('Rise Time');
figure('Name','Montecarlo: Delta lateral','NumberTitle','off'), histogram(delta_lat_mc,100), grid on, title('Delta lateral');
xlabel('Delta lateral');














