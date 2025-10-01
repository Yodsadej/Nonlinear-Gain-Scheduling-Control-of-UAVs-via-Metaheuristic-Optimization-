% PID Tuning 
%-------->  Va
% |
% |
% |
% V
% altitude

function [fp,f,g,fout1,fout2,gr,SysLIC0,SysLIC1,SysLIC]=UAVTuningResult(x,sysL,ref,X0)
warning('off')
% load('RefVelocity0.mat')
% sysL=sys;
try
%% construct the interconnection
% simulate actuator and its rate
ag=25;
Act=[tf([0 ag],[1 ag])
    tf([ag 0],[1 ag])];

Act.u='u';
Act.y={'dE','dE_rate'};

% Va2Gamma=pid(-0.1,-0.1,0);
% Gamma2dE=pid(-17.2161065252514,-8.17127344979209,-8.96022389857252,1/138.399392031606);


Va2Gamma=pid(x(1),x(2),x(3),x(4));
Gamma2dE=pid(x(5),x(6),x(7),x(8));

sysL.u={'dE','dTh'};
sysL.y={'Va','gamma','alpha','q'};

Va2Gamma.u='ErrorVa'; %command speed
Va2Gamma.y='Cgamma'; % command gamma

Gamma2dE.u='Errorgamma'; 
Gamma2dE.y='u'; % to actuator


sum1=sumblk('ErrorVa=CVa-Va');
sum2=sumblk('Errorgamma=Cgamma-gamma');

opt = connectOptions('Simplify',false);
% system close loop interconnection

SysLIC0=connect(sysL,Act,Gamma2dE,...
        {'Errorgamma'},{'Va','gamma','alpha','q','dE','dE_rate'},opt);

SysLIC1=connect(sysL,Act,Gamma2dE,Va2Gamma,sum2,...
        {'ErrorVa'},{'Va','gamma','alpha','q','dE','dE_rate'},opt);

SysLIC=connect(sysL,Act,Va2Gamma,Gamma2dE,sum1,sum2,...
        {'CVa'},{'Va','gamma','alpha','q','dE','dE_rate','Cgamma'},opt);
        % {'CVa','dTh'},{'Va','gamma','alpha','q','dE','dE_rate'},opt);


% sum1=sumblk('yG(1)=yd(1)+theta');

% step(SysLIC(:,1,1,1))
% save RefVelocity ref t

% L=Gamma2dE*sysL(2,1,1,1);
% L2=sysL(2,1,1,1)*Gamma2dE;
% diskmargin(Gamma2dE*sysL(2,1,1,1)*Act(1,1))
% % diskmargin(SysLIC0(2,1,1,1))
% Gamma2dE*sysL(1,1,1,1)

%% cheak stability
if isstable(SysLIC)==0
    f=1e15;
    fp=1e15;
    g=1e15*ones(3,1);
    return;
end

%% step evaluation

finalvalue=max(ref);
% opt = stepDataOptions('StepAmplitude',finalvalue);
opt = RespConfig('InitialState',[X0;0;0;0],'StepAmplitude',finalvalue);
[y,t]=step(SysLIC(:,1),100,opt); %y(1),y(2)
% [y,~]=lsim(SysLIC(:,1,1,1),ref0,t,X0);
info=stepinfo(y(:,1),t,finalvalue); % Va
ctrl_eff=sqrt(mean((y(:,5).^2)));
% if err

%% set constaint
elevator_limit=10; %+-10
elevator_rate_limit=50;
gamma_command_limit=0.2; %0.2 rad, around 11.46 degree
f10=info.Overshoot;
% g1=(info.Overshoot-5)/5;
g1=(max(abs(y(:,5)))-elevator_limit)/elevator_limit; %% aeileron defection limit
g2=(max(abs(y(:,6)))-elevator_rate_limit)/elevator_rate_limit; % y6 already is a rate
g3=(max(abs(y(:,7)))-gamma_command_limit)/gamma_command_limit; % 

gr1=max(abs(y(:,5)));
gr2=max(abs(y(:,6)));
gr3=max(abs(y(:,7)));
% g3=(max(abs(diff(y(:,6)))./t(2))-elevator_rate_limit)/elevator_rate_limit;% control effort rate limit

%% set obj
f1=info.SettlingTime;
% f2=sqrt(mean((y(0.9*end:end,1)-finalvalue).^2)); %% sse
f2=sqrt(mean((y(:,1)-finalvalue).^2)); %% sse
f3=ctrl_eff; %% control effort

%% evaluate margin
% openloop gamma to elevator
S0=allmargin(SysLIC0(2,1));
[DM0,~]=diskmargin(SysLIC0(2,1));

% closeloop gamma to elevator--> openloop Velocity to gamma
S1=allmargin(SysLIC1(1,1));
[DM1,~]=diskmargin(SysLIC1(1,1));

f4=1/DM0.DiskMargin;
f5=1/abs(min(S0.GainMargin));
f6=1/abs(min(S0.PhaseMargin));
f7=1/DM1.DiskMargin;
f8=1/abs(min(S1.GainMargin));
f9=1/abs(min(S1.PhaseMargin));

fr4=DM0.DiskMargin;
fr5=abs(min(S0.GainMargin));
fr6=abs(min(S0.PhaseMargin));
fr7=DM1.DiskMargin;
fr8=abs(min(S1.GainMargin));
fr9=abs(min(S1.PhaseMargin));

catch me    
    f=1e15;
    fp=1e15;
    g=1e15*ones(3,1);% 5 constains
%     load logdata.mat
%     err=err+1;
%     logdata{err}=me;
%     save logdata.mat err logdata;
     return;
end


%% cleack Nan
if isnan(g1)
    g1=1e10;
end
if isnan(g2)
    g2=1e10;
end
if isnan(g3)
    g3=1e10;
end
if isnan(f1)
    f1=1e10;
end
if isnan(f2)
    f2=1e10;
end
if isnan(f3)
    f3=1e10;
end

g=[g1;g2;g3];
f=0.1*f1+10*f2+f3+f4+f5+f6+f7+f8+f9;
fout1=[f1 f2 f3 fr4 fr5 fr6 fr7 fr8 fr9]; % result from objective
fout2=[f1 f2 f3 f4 f5 f6 f7 f8 f9]; % Margin
gr=[g1;g2;g3];
% fp=f+g;
fp=fpenal(f,g);

function fp=fpenal(f,g)
fp=f+100*max(0,max(g));







