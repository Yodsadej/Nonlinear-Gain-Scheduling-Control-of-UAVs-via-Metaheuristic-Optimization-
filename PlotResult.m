clear all;close all;clc;
load SimulationData.mat
load Usampling.mat
out.Va=out.Va(:,:)';
nt=length(out.Va);
% Vc=ones(nt,1)*out.Vc;


%% plot altitude and velocity
figure()
s1=subplot(2,1,1)
p=plot(out.tout,out.hc,out.tout,out.h);

set(gcf,'Position',[0,0,800,500])

ylabel('Altitude (m)');


grid on
p(1).LineWidth=1.5;
p(1).LineStyle='--';
p(2).LineWidth=1.25;

legend('h_{cmd}','h')
ylim([0,9000])
xticklabels("")


s2=subplot(2,1,2)
p2=plot(out.tout,out.Vc,out.tout,out.Va(:,:))
ylabel('Velocity (m/s)');
xlabel('time (sec)');

p2(1).LineWidth=1.5;
p2(1).LineStyle='--';
p2(2).LineWidth=1.25;
grid on

legend('V_{cmd}','V_a')
legend('Location','southeast')

pos1 = get(s1, 'Position') % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0]
set(s2, 'Position',new_pos1 ) % set new position of current sub - plot

%% plot flight status mode and throttle command end elevator command and fligh path angle


figure()
s3=subplot(2,1,1)
p3=plot(out.tout,out.dE);

set(gcf,'Position',[0,0,800,500])

ylabel('Altitude (m)');


grid on
p3(1).LineWidth=1.5;

legend('dE')
% ylim([0,9000])
xticklabels("")


s4=subplot(2,1,2)
p4=plot(out.tout,out.dT);
ylabel('Velocity (m/s)');
xlabel('time (sec)');

p4(1).LineWidth=1.5;
grid on

legend('V_{cmd}','V_a')
legend('Location','southeast')

pos1 = get(s1, 'Position') % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0]
set(s2, 'Position',new_pos1 ) % set new position of current sub - plot


%% combine figure "experimental"
close all
gap=-0.115;
figure()
s1=subplot(6,1,1)
for j=1:50
    eval(['plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.Va,''Color'', [0.7 0.7 0.7])'])
    hold on
end
p=plot(out.tout,out.Va(:,:),out.tout,out.Vc);
xlim([0,4750])
ylim([35,65])
set(gcf,'Position',[0,0,1200,1500])
ylabel('Velocity (m/s)');
legend(p,'V_a','V_{cmd}')
% legend('Location','southeast')
fontsize(14,"points")

grid on
p(1).LineWidth=1.5;
p(1).Color='[0 0.4470 0.7410]';
% p(1).LineStyle='--';
p(2).LineWidth=1.75;
p(2).LineStyle='--';
p(2).Color='red';
xticklabels("")


s2=subplot(6,1,2)
for j=1:50
    eval(['plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.h,''Color'', [0.7 0.7 0.7])'])
    hold on
end
p2=plot(out.tout,out.h,out.tout,out.hc);
xlim([0,4750])
ylabel('Altitude (m)');
legend(p2,'h','h_{cmd}')

p2(1).LineWidth=1.5;
p2(1).Color='[0 0.4470 0.7410]';
p2(2).LineWidth=1.75;
p2(2).LineStyle='--';
p2(2).Color='red';
xticklabels("")

ylim([0,9000])
grid on
xticklabels("")
fontsize(14,"points")

pos1 = get(s1, 'Position') % gives the position of current sub-plot
new_pos1 = pos1 +[0 gap 0 0]
set(s2, 'Position',new_pos1 ) % set new position of current sub - plot


%%%%%%%%%%%%%%%%%%%%%%%%%%

s3=subplot(6,1,3)
for j=1:50
    eval(['plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.gamma,''Color'', [0.7 0.7 0.7])'])
    hold on
end
p3=plot(out.tout,out.gamma,out.tout,out.GammaCommand);
xlim([0,4750])
ylabel('\gamma (deg)');
legend(p3,'\gamma','\gamma_{cmd}')

fontsize(14,"points")
% ylim([-2,10])

p3(1).LineWidth=1.5;
p3(1).Color='[0 0.4470 0.7410]';
p3(2).LineWidth=1.75;
p3(2).LineStyle='--';
p3(2).Color='red';
xticklabels("")
grid on
xticklabels("")

pos2 = get(s2, 'Position') % gives the position of current sub-plot
new_pos3 = pos2 +[0 gap 0 0]
set(s3, 'Position',new_pos3 ) % set new position of current sub - plot

%%%%%%%%%%%%%%%%%%%%%%%%%%
s4=subplot(6,1,4)
for j=1:50
    eval(['plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.dE,''Color'', [0.7 0.7 0.7])'])
    hold on
end
p4=plot(out.tout,out.dE);
xlim([0,4750])
ylabel('Elevator Deflection (deg)');
fontsize(14,"points")
ylim([-2,10])

p4(1).LineWidth=1.5;
p4(1).Color='[0 0.4470 0.7410]';
grid on
xticklabels("")


pos3 = get(s3, 'Position') % gives the position of current sub-plot
new_pos4 = pos3 +[0 gap 0 0]
set(s4, 'Position',new_pos4 ) % set new position of current sub - plot

%%%%%%%%%%%%%%%%%%%%%%%%%%
s5=subplot(6,1,5)
for j=1:50
    eval(['plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.dT,''Color'', [0.7 0.7 0.7])'])
    hold on
end
% p4=plot(out.tout,out.dT,out.tout,0.5.*out.status)
p5=plot(out.tout,out.dT);
xlim([0,4750])

p5(1).LineWidth=1.5;
p5(1).Color='[0 0.4470 0.7410]';
% p4(2).LineWidth=1.5;

ylabel('Throttle (%)');
fontsize(14,"points")
ylim([-0.1,1.1])
grid on
xticklabels("")

pos4 = get(s4, 'Position') % gives the position of current sub-plot
new_pos5 = pos4 +[0 gap 0 0]
set(s5, 'Position',new_pos5 ) % set new position of current sub - plot

%%%%%%%%%%%%%%%%%%%%%%%%%%
s6=subplot(6,1,6)
for j=1:50
    eval(['plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.status,''Color'', [0.7 0.7 0.7])'])
    hold on
end
p6=plot(out.tout,out.status);
xlim([0,4750])

p6(1).LineWidth=2.5;
p6(1).Color='[0 0.4470 0.7410]';
grid on

pos5 = get(s5, 'Position') % gives the position of current sub-plot
new_pos6 = pos5 +[0 gap 0 0]
set(s6, 'Position',new_pos6 ) % set new position of current sub - plot

ylabel('Flight Mode');
xlabel('time (sec)');
fontsize(14,"points")




%%
figure

% p=plot(out.tout,Vc,out.tout,out.Va(:,:))
% p2=plot(out.tout,out.hc,out.tout,out.h);
% p3=plot(out.tout,out.dE);
% p4=plot(out.tout,out.dT);

t1start=1450;
t1end=1550;
t2start=1900;
t2end=2600;

t3start=3300;
t3end=3900;
% t4start=4400;
% t4end=4500;

t1=out.tout(out.tout>t1start & out.tout<t1end);
Vc_cut1=out.Vc(out.tout>t1start & out.tout<t1end);
Va_cut1=out.Va(out.tout>t1start & out.tout<t1end);
gammac_cut1=out.GammaCommand(out.tout>t1start & out.tout<t1end);
gamma_cut1=out.gamma(out.tout>t1start & out.tout<t1end);
hc_cut1=out.hc(out.tout>t1start & out.tout<t1end);
h_cut1=out.h(out.tout>t1start & out.tout<t1end);
de_cut1=out.dE(out.tout>t1start & out.tout<t1end);
dT_cut1=out.dT(out.tout>t1start & out.tout<t1end);


t2=out.tout(out.tout>t2start & out.tout<t2end);
Vc_cut2=out.Vc(out.tout>t2start & out.tout<t2end);
Va_cut2=out.Va(out.tout>t2start & out.tout<t2end);
gammac_cut2=out.GammaCommand(out.tout>t2start & out.tout<t2end);
gamma_cut2=out.gamma(out.tout>t2start & out.tout<t2end);
hc_cut2=out.hc(out.tout>t2start & out.tout<t2end);
h_cut2=out.h(out.tout>t2start & out.tout<t2end);
de_cut2=out.dE(out.tout>t2start & out.tout<t2end);
dT_cut2=out.dT(out.tout>t2start & out.tout<t2end);


t3=out.tout(out.tout>t3start & out.tout<t3end);
Vc_cut3=out.Vc(out.tout>t3start & out.tout<t3end);
Va_cut3=out.Va(out.tout>t3start & out.tout<t3end);
gammac_cut3=out.GammaCommand(out.tout>t3start & out.tout<t3end);
gamma_cut3=out.gamma(out.tout>t3start & out.tout<t3end);
hc_cut3=out.hc(out.tout>t3start & out.tout<t3end);
h_cut3=out.h(out.tout>t3start & out.tout<t3end);
de_cut3=out.dE(out.tout>t3start & out.tout<t3end);
dT_cut3=out.dT(out.tout>t3start & out.tout<t3end);

% t4=out.tout(out.tout>t4start & out.tout<t4end);
% Vc_cut4=out.Vc(out.tout>t4start & out.tout<t4end);
% Va_cut4=out.Va(out.tout>t4start & out.tout<t4end);
% hc_cut4=out.hc(out.tout>t4start & out.tout<t4end);
% h_cut4=out.h(out.tout>t4start & out.tout<t4end);
% de_cut4=out.dE(out.tout>t4start & out.tout<t4end);
% dT_cut4=out.dT(out.tout>t4start & out.tout<t4end);


gapcx=0.23;
gapcy=-0.18;

n_section=3;



for i=1:n_section
    eval(['sc' num2str(i) '=subplot(4,' num2str(n_section) ',' num2str(i) ');'])
    hold on
    for j=1:50
        eval(['cc=plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.Va,''Color'', [0.7 0.7 0.7]);'])
        hold on
    end
    eval(['c' num2str(i) '=plot(t' num2str(i) ',Va_cut' num2str(i) ',t' num2str(i) ',Vc_cut' num2str(i) ');'])
    eval(['c' num2str(i) '(1).LineWidth=1.5;'])
    eval(['c' num2str(i) '(1).Color=''[0 0.4470 0.7410]'';'])
    eval(['c' num2str(i) '(2).LineStyle=''--'';'])
    eval(['c' num2str(i) '(2).LineWidth=1.25;'])
    eval(['c' num2str(i) '(2).Color=''red'';'])

    grid on
    ylim([35,65])
    
    if i==1
        set(gcf,'InnerPosition',[0,0,1500,1000])
        ylabel('Velocity (m/s)','fontsize',12);
    end

    if i~=1
        eval(['posc' num2str(i-1) ' = get(sc' num2str(i-1) ', ''Position'');'])
        eval(['new_posc' num2str(i) ' = posc' num2str(i-1) ' +[gapcx 0 0 0];'])
        eval(['set(sc' num2str(i) ', ''Position'',new_posc' num2str(i) ');'])
        yticklabels("")
        
    end

    eval(['xlim([t' num2str(i) 'start,t' num2str(i) 'end]);'])
    xticklabels("")
end
legend(c3,'V_a','V_{cmd}','fontsize',12)
       % legend(cc2,'V_a_{Uncertainty}','fontsize',12)
       % legend('Location','southwest')


posc1 = get(sc1, 'Position');
new_posc4 = posc1 +[0 gapcy 0 0];

for i=n_section+1:2*n_section
    eval(['sc' num2str(i) '=subplot(4,' num2str(n_section) ',' num2str(i) ');'])
    hold on
    for j=1:50
        eval(['cc=plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.gamma,''Color'', [0.7 0.7 0.7]);'])
        hold on
    end
    eval(['c' num2str(i) '=plot(t' num2str(i-n_section) ',gamma_cut' num2str(i-n_section) ',t' num2str(i-n_section) ',gammac_cut' num2str(i-n_section) ');'])
    eval(['c' num2str(i) '(1).LineWidth=1.5;'])
    eval(['c' num2str(i) '(1).Color=''[0 0.4470 0.7410]'';'])
    eval(['c' num2str(i) '(2).LineStyle=''--'';'])
    eval(['c' num2str(i) '(2).LineWidth=1.25;'])
    eval(['c' num2str(i) '(2).Color=''red'';'])
    grid on
    ylim([-25,25])
    
    if i==4
        % xlim=([2240,2310]);
        ylabel('\gamma (deg)','fontsize',12);
        set(sc4, 'Position',new_posc4 )
    end

    if i~=4
        eval(['posc' num2str(i-1) ' = get(sc' num2str(i-1) ', ''Position'');'])
        eval(['new_posc' num2str(i) ' = posc' num2str(i-1) ' +[gapcx 0 0 0];'])
        eval(['set(sc' num2str(i) ', ''Position'',new_posc' num2str(i) ');'])
        yticklabels("")
    end

    % if i==2*n_section
    %     legend('\gamma_{cmd}','\gamma','fontsize',12)
    % end

    eval(['xlim([t' num2str(i-n_section) 'start,t' num2str(i-n_section) 'end]);'])
    xticklabels("")
end
legend(c6,'\gamma','\gamma_{cmd}','fontsize',12)
% legend([c6;cc]','\gamma','\gamma_{cmd}','\gamma_{Usampling}','fontsize',12)


posc4 = get(sc4, 'Position');
new_posc7 = posc4 +[0 gapcy 0 0];


for i=2*n_section+1:3*n_section
    eval(['sc' num2str(i) '=subplot(4,' num2str(n_section) ',' num2str(i) ');'])
    hold on
    for j=1:50
        eval(['cc=plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.dE,''Color'', [0.7 0.7 0.7]);'])
        hold on
    end
    eval(['c' num2str(i) '=plot(t' num2str(i-2*n_section) ',de_cut' num2str(i-2*n_section) ');'])
    eval(['c' num2str(i) '(1).LineWidth=1.5;'])
    eval(['c' num2str(i) '(1).Color=''[0 0.4470 0.7410]'';'])
    ylim([-2,10])
    grid on

if i==7
    % xlim=([2240,2310]);
    ylabel('Elevator Deflection (deg)','fontsize',12);
    set(sc7, 'Position',new_posc7 )
end
if i~=7
    eval(['posc' num2str(i-1) ' = get(sc' num2str(i-1) ', ''Position'');'])
    eval(['new_posc' num2str(i) ' = posc' num2str(i-1) ' +[gapcx 0 0 0];'])
    eval(['set(sc' num2str(i) ', ''Position'',new_posc' num2str(i) ');'])
    yticklabels("")
end
eval(['xlim([t' num2str(i-2*n_section) 'start,t' num2str(i-2*n_section) 'end]);'])
xticklabels("")
% xtickangle(90)
% eval(['xlabel(''Time interval ' num2str(i-2*n_section) ' (sec)'',''fontsize'',12)'])
end

posc7 = get(sc7, 'Position');
new_posc10 = posc7 +[0 gapcy 0 0];

for i=3*n_section+1:4*n_section
    eval(['sc' num2str(i) '=subplot(4,' num2str(n_section) ',' num2str(i) ');'])
    hold on
    for j=1:50
        eval(['plot(Usamplimg' num2str(j) '.tout,Usamplimg' num2str(j) '.dT,''Color'', [0.7 0.7 0.7])'])
        hold on
    end
    eval(['c' num2str(i) '=plot(t' num2str(i-3*n_section) ',dT_cut' num2str(i-3*n_section) ');'])
    eval(['c' num2str(i) '(1).LineWidth=1.5;'])
    eval(['c' num2str(i) '(1).Color=''[0 0.4470 0.7410]'';'])
    ylim([-0.1,1.1])
    grid on

if i==10
    % xlim=([2240,2310]);
    ylabel('Throttle (%)','fontsize',12);
    set(sc10, 'Position',new_posc10 )
end
if i~=10
    eval(['posc' num2str(i-1) ' = get(sc' num2str(i-1) ', ''Position'');'])
    eval(['new_posc' num2str(i) ' = posc' num2str(i-1) ' +[gapcx 0 0 0];'])
    eval(['set(sc' num2str(i) ', ''Position'',new_posc' num2str(i) ');'])
    yticklabels("")
end
eval(['xlim([t' num2str(i-3*n_section) 'start,t' num2str(i-3*n_section) 'end]);'])
xtickangle(90)
eval(['xlabel(''Time interval ' num2str(i-3*n_section) ' (sec)'',''fontsize'',12)'])

end


%% compare with and without descending command


load WithdescendingCtrl
load NodescendingCtrl
t_comp_start=2000;


tcomp1=WithdescendingCtrl.tout(WithdescendingCtrl.tout>t_comp_start);
Vc_comp1=WithdescendingCtrl.Vc(WithdescendingCtrl.tout>t_comp_start);
Va_comp1=WithdescendingCtrl.Va(WithdescendingCtrl.tout>t_comp_start);
% gammac_cut1=out.GammaCommand(out.tout>t1start & out.tout<t1end);
% gamma_cut1=out.gamma(out.tout>t1start & out.tout<t1end);
hc_comp1=WithdescendingCtrl.hc(WithdescendingCtrl.tout>t_comp_start);
h_comp1=WithdescendingCtrl.h(WithdescendingCtrl.tout>t_comp_start);


tcomp2=NodescendingCtrl.tout(NodescendingCtrl.tout>t_comp_start);
Vc_comp2=NodescendingCtrl.Vc(NodescendingCtrl.tout>t_comp_start);
Va_comp2=NodescendingCtrl.Va(NodescendingCtrl.tout>t_comp_start);
% gammac_cut1=out.GammaCommand(out.tout>t1start & out.tout<t1end);
% gamma_cut1=out.gamma(out.tout>t1start & out.tout<t1end);
hc_comp2=NodescendingCtrl.hc(NodescendingCtrl.tout>t_comp_start);
h_comp2=NodescendingCtrl.h(NodescendingCtrl.tout>t_comp_start);


figure()
s1=subplot(2,1,1)
p=plot(tcomp1,h_comp1,tcomp2,h_comp2,tcomp1,hc_comp1);

set(gcf,'Position',[0,0,800,500])

ylabel('Altitude (m)');

grid on

p(1).LineWidth=1.5;
p(2).LineWidth=1.5;
% p(1).LineStyle='--';
p(3).LineWidth=1.75;
p(3).LineStyle='--';
p(3).Color='red';
xticklabels("")

legend('dual controller','single controller','h_{cmd}')
ylim([0,9000])
xlim([t_comp_start,4750])
xticklabels("")

s2=subplot(2,1,2)
p2=plot(tcomp1,Va_comp1,tcomp2,Va_comp2,tcomp1,Vc_comp1);
ylabel('Velocity (m/s)');
xlabel('time (sec)');
xlim([t_comp_start,4750])

p2(1).LineWidth=1.5;
p2(2).LineWidth=1.5;
p2(3).LineWidth=1.75;
p2(3).LineStyle='--';
p2(3).Color='red';

grid on

legend('dual controller','single controller','V_{cmd}')


pos1 = get(s1, 'Position') % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0]
set(s2, 'Position',new_pos1 ) % set new position of current sub - plot

%% uncertainty sampling
figure()
s1=subplot(2,1,1)
for i=1:50
    eval(['plot(Usamplimg' num2str(i) '.tout,Usamplimg' num2str(i) '.h,''Color'', [0.5 0.5 0.5])'])
    hold on
end

p=plot(WithdescendingCtrl.tout,WithdescendingCtrl.h,WithdescendingCtrl.tout,WithdescendingCtrl.hc,'Color',  [0 0.4470 0.7410]);
set(gcf,'Position',[0,0,800,500])

ylabel('Altitude (m)');

grid on
p(1).LineWidth=1.5;
% p(2).LineWidth=1.5;
% p(1).LineStyle='--';
p(2).LineWidth=1.75;
p(2).LineStyle='--';
p(2).Color='red';
xticklabels("")

% legend('h_{cmd}','h_w')
ylim([0,9000])
xlim([0,4750])
xticklabels("")
legend('h_{cmd}','h_w')

s2=subplot(2,1,2)

for i=1:50
    eval(['plot(Usamplimg' num2str(i) '.tout,Usamplimg' num2str(i) '.Va,''Color'', [0.5 0.5 0.5])'])
    hold on
end

p2=plot(WithdescendingCtrl.tout,WithdescendingCtrl.Vc,WithdescendingCtrl.tout,WithdescendingCtrl.Va,'Color',  [0 0.4470 0.7410]);
ylabel('Velocity (m/s)');
xlabel('time (sec)');
ylim([35,65])
xlim([0,4750])

p2(2).LineWidth=1.5;
% p2(2).LineWidth=1.5;
p2(1).LineWidth=1.75;
p2(1).LineStyle='--';
p2(1).Color='red';

grid on

legend('V_{cmd}','V_a')

pos1 = get(s1, 'Position') % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0]
set(s2, 'Position',new_pos1 ) % set new position of current sub - plot
% load('Usampling.mat')

% load WithdescendingCtrl
% plot(WithdescendingCtrl.tout,WithdescendingCtrl.Vc,'r--')
% hold on
% for i=1:50
%     eval(['plot(Usamplimg' num2str(i) '.tout,Usamplimg' num2str(i) '.Va,''Color'', [0.5 0.5 0.5])'])
%     hold on
% end




% WithdescendingCtrl=out;

% NodescendingCtrl=out;
% save NodescendingCtrl NodescendingCtrl
% save WithdescendingCtrl WithdescendingCtrl


