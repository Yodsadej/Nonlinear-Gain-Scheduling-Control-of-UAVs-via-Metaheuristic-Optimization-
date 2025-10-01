close all;clear all;clc;

load('EUAV_ss05.mat')
load('EUAV_ssL5.mat')
load('RefVelocity.mat')

nrun=2;

% --------> velocity
% |
% |
% |
% V
% Altitude

% foutput=[char(objfun(k,:)) '_' char(algo(i,:)) ' run ' num2str(j) '_' num2str(h) '_' num2str(V)];
for j0=1:size(sysL,3) % altitude 
    for j1=1:size(sysL,4) % speed
        sys=sysL(:,:,j0,j1);
        X0=[Xtt(:,j1,:,j0);0;0];
        disp(X0(1))
        disp(sysL(:,:,j0,j1).SamplingGrid)
        h=sysL(:,:,j0,j1).SamplingGrid.h;
        V=sysL(:,:,j0,j1).SamplingGrid.V;

        ref0=ones(length(t),1);
        ref=(V-10)*ref0;
        xmint=[];fpmint=[];gmint=[];
        for i=1:nrun
            % ['load(''UAVTuning_SOLSHADE run ' num2str(i) '_' num2str(h) '_' num2str(V) '.mat;'')']
            eval(['load(''UAVTuning_SOLSHADE run ' num2str(i) '_' num2str(h) '_' num2str(V) '.mat'')'])
            xmint=[xmint;xmin];
            fpmint=[fpmint;fpmin];
            gmint=[gmint,gmin];
        end
        [minValue,index]=min(fpmint);
        gmintVH(:,:,j0,j1)=gmint(:,index);
        xminVH(j0,j1,:)=xmint(index,:);
        fpminVH(j0,j1)=minValue;
        Utrim(j0,j1,:)=Triminfot(:,j0,1,j1)';
        
        % [~,~,~,fout1,fout2,gr]=UAVTuningUncertaintySampling(xmint(index,:),sys,ref,X0);
        [~,~,~,fout1,fout2,gr,SysLIC0(:,:,j0,j1),SysLIC1(:,:,j0,j1),SysLIC(:,:,j0,j1)]=UAVTuningResult(xmint(index,:),sys,ref,X0);
        

        ST(j0,j1)=fout1(1); %settling time
        SSE(j0,j1)=fout1(2); %settling time
        Ctrl_eff(j0,j1)=fout1(3); %settling time
        DM0(j0,j1)=fout1(4); % innerloop disk margin
        GainMargin0(j0,j1)=fout1(5);  % innerloop gain margin
        PhaseMargin0(j0,j1)=fout1(6); % innerloop phase margin
        DM2(j0,j1)=fout1(7); % outerloop disk margin
        GainMargin2(j0,j1)=fout1(8);  % outerloop gain margin
        PhaseMargin2(j0,j1)=fout1(9); % outerloop phase margin
    end
end

%% plot data
Vsampling=30:5:80;
Altitude=0:1000:10000;
Altitude=num2cell(Altitude);
% legendEntries = {'dev1' 'dev2' 'dev3' 'dev4' 'dev5' 'dev6'};
Speeds=num2cell(Vsampling);


%% plot data surface
% facecolor=[0.6,0.6,0.6]; % grey
% facecolor=[255/256, 192/256, 203/256]; % purple
% facecolor=[233, 30, 99 ]./256; % pink
facecolor=[1, 0, 0 ];
x0=10;
y0=10;
width=550;
height=550;

% facecolor=[1,1,1]; % 

% subplot(1,3,1)
figure()
surf(ST)
view(-45,45)
% pbaspect([1 1 1])
xlim([1,11]);
ylim([1,11]);
set(gca,'XTickLabel',Speeds,'YTickLabel',Altitude,'Ytick',1:11,'XTick',1:11)
hold on
% plot(mean(ST),'-x','Color','r')
colorbar("north")
grid on
xlabel('Velocity')
ylabel('Altitude')
zlabel('second')

set(gcf,'position',[x0,y0,width,height])
% title('a) Settling time','FontSize',12) 

% subplot(1,3,2)
figure
surf(SSE)
view(-45,45)
hold on
% plot(mean(SSE),'-x','Color','r')
xlim([1,11]);
ylim([1,11]);
set(gca,'XTickLabel',Speeds,'YTickLabel',Altitude,'Ytick',1:11,'XTick',1:11)
colorbar("north")

% set(gca,'XTickLabel',Speeds)
grid on
xlabel('Velocity')
ylabel('Altitude')
zlabel('RMS')
% title('b) Steady state error','FontSize',12) 
set(gcf,'position',[x0,y0,width,height])

% subplot(1,3,3)
figure()
surf(Ctrl_eff)
view(-45,45)
hold on
xlim([1,11]);
ylim([1,11]);
% plot(mean(Ctrl_eff),'-x','Color','r')
set(gca,'XTickLabel',Speeds,'YTickLabel',Altitude,'Ytick',1:11,'XTick',1:11)
colorbar("north")
grid on

xlabel('Velocity')
ylabel('Altitude')
zlabel('RMS')
% title('c) Control effort','FontSize',12) 
set(gcf,'position',[x0,y0,width,height])


%% Bode plot

w = logspace(-4,4,1000);
for i=1:11
    for j=1:11
        [mag,phase,wout] =bode(SysLIC1(1,1,i,j),w); % loop transfer
        % [mag2,phase2,wout2] =bode(SysLIC(1,1,i,j),w);
        s1=subplot(2,1,1);
        semilogx(wout, 20*log10(squeeze(mag)), 'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        hold on
        % semilogx(wout2, 20*log10(squeeze(mag2)),'--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Magnitude (dB)');
        xticklabels("")
        
        grid
        s2=subplot(2,1,2);
        semilogx(wout, squeeze(phase), 'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        ylabel('Phase (deg)');
        hold on
        grid
    end
end

for i=1:11
    for j=1:11
        [mag3,phase3,wout3] =bode(1/(1+SysLIC1(1,1,i,j)),w); % sensivity
        % [mag2,phase2,wout2] =bode(SysLIC(1,1,i,j),w);
        s1=subplot(2,1,1);
        semilogx(wout3, 20*log10(squeeze(mag3)), 'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        hold on
        % semilogx(wout2, 20*log10(squeeze(mag2)),'--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Magnitude (dB)');
        xticklabels("")
        
        grid
        s2=subplot(2,1,2);
        semilogx(wout3, squeeze(phase3), 'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        ylabel('Phase (deg)');
        hold on
        grid
    end
end

for i=1:11
    for j=1:11
        [mag2,phase2,wout2] =bode(SysLIC(1,1,i,j),w);  % transfer function
        subplot(2,1,1)
        semilogx(wout2, 20*log10(squeeze(mag2)),'--k', 'LineWidth',1.25)
        xticklabels("")

        subplot(2,1,2)
        semilogx(wout2, squeeze(phase2),'--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        xlabel('Frequency (rad/s)');
        bandwidth_loop2(i,j)=bandwidth(SysLIC(1,1,i,j));
    end
end
% legend('Sensivity function','Comprementary')


pos1 = get(s1, 'Position'); % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0];
set(s2, 'Position',new_pos1 );% set new position of current sub - plot

min(min(bandwidth_loop2))
max(max(bandwidth_loop2))
mean(mean(bandwidth_loop2))
%%

w = logspace(-4,4,1000);
figure()
for i=1:11
    for j=1:11
        [mag0,phase0,wout0] =bode(SysLIC0(2,1,i,j),w); % loop transfer
        % [mag2,phase2,wout2] =bode(SysLIC(1,1,i,j),w);
        s1=subplot(2,1,1);
        semilogx(wout0, 20*log10(squeeze(mag0)), 'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        hold on
        % semilogx(wout2, 20*log10(squeeze(mag2)),'--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Magnitude (dB)');
        xticklabels("")
        
        grid
        s2=subplot(2,1,2);
        semilogx(wout0, squeeze(phase0), 'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Phase (deg)');
        xlabel('Frequency (rad/s)');
        hold on
        grid
    end
end

for i=1:11 
    for j=1:11
        [mag02,phase02,wout02] =bode(feedback(SysLIC0(2,1,i,j),1),w); % transfer function
        % [mag2,phase2,wout2] =bode(SysLIC(1,1,i,j),w);
        s1=subplot(2,1,1);
        semilogx(wout02, 20*log10(squeeze(mag02)), '--k', 'LineWidth',1.25)
        hold on
        % semilogx(wout2, 20*log10(squeeze(mag2)),'--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Magnitude (dB)');
        xticklabels("")
        
        grid
        s2=subplot(2,1,2);
        semilogx(wout02, squeeze(phase02), '--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Phase (deg)');
        xlabel('Frequency (rad/s)');
        hold on
        grid

        bandwidth_loop1(i,j)=bandwidth(feedback(SysLIC0(2,1,i,j),1));
    end
end

for i=1:11 
    for j=1:11
        [mag03,phase03,wout03] =bode(1/(1+SysLIC0(2,1,i,j)),w); % sensivity function
        % [mag2,phase2,wout2] =bode(SysLIC(1,1,i,j),w);
        s1=subplot(2,1,1);
        semilogx(wout03, 20*log10(squeeze(mag03)), 'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        hold on
        % semilogx(wout2, 20*log10(squeeze(mag2)),'--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Magnitude (dB)');
        xticklabels("")
        
        grid
        s2=subplot(2,1,2);
        semilogx(wout03, squeeze(phase03), 'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Phase (deg)');
        xlabel('Frequency (rad/s)');
        hold on
        grid
    end
end


pos1 = get(s1, 'Position'); % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0];
set(s2, 'Position',new_pos1 );% set new position of current sub - plot

min(min(bandwidth_loop1))
max(max(bandwidth_loop1))
mean(mean(bandwidth_loop1))

%% for flightpath tracking
for j0=1:size(sysL,3)
    for j1=1:size(sysL,4)
        sys=sysL(:,:,j0,j1);
        h=sysL(:,:,j0,j1).SamplingGrid.h;
        V=sysL(:,:,j0,j1).SamplingGrid.V;
        xmint=[];fpmint=[];
        for i=1:nrun
            % ['load(''UAVTuningGamma_SOLSHADE run ' num2str(i) '_' num2str(h) '_' num2str(V) '.mat;'')']
            eval(['load(''UAVTuningGamma_SOLSHADE run ' num2str(i) '_' num2str(h) '_' num2str(V) '.mat'')'])
            xmint=[xmint;xmin];
            fpmint=[fpmint;fpmin];
        end
        [minValue,index]=min(fpmint);
        xminVH_gamma(j0,j1,:)=xmint(index,1:4);
        fpminVH_gamma(j0,j1)=minValue;
        % Utrim(j0,j1,:)=Triminfot(:,j0,1,j1)';

        [~,~,g,fout1,fout2,SysLICg0(:,:,j0,j1),SysLICg(:,:,j0,j1)]=UAVTuningGammaResult(xmint(index,:),sys);
        ST(j0,j1)=fout1(1); % settling time
        SSE(j0,j1)=fout1(2); % steady state error
        Ctrl_eff(j0,j1)=fout1(3); % control effort
        DM0(j0,j1)=fout1(4); % innerloop disk margin
        GainMargin0(j0,j1)=fout1(5);  % innerloop gain margin
        PhaseMargin0(j0,j1)=fout1(6); % innerloop phase margin
    end
end

d2r=pi/180;
pid_h2gam=pid(0.2,0.005);lowpass=tf(1,[2,1]);

CL=feedback(d2r*pid_h2gam*lowpass*SysLICg(2,1,:,:),1);


%% plot data surface
% facecolor=[0.6,0.6,0.6]; % grey
% facecolor=[255/256, 192/256, 203/256]; % purple
% facecolor=[233, 30, 99 ]./256; % pink
facecolor=[1, 0, 0 ];
x0=10;
y0=10;
width=550;
height=550;

% facecolor=[1,1,1]; % 

% subplot(1,3,1)
figure()
surf(ST)
view(45,45)
% pbaspect([1 1 1])
xlim([1,11]);
ylim([1,11]);
set(gca,'XTickLabel',Speeds,'YTickLabel',Altitude,'Ytick',1:11,'XTick',1:11)
hold on
% plot(mean(ST),'-x','Color','r')
colorbar("north")
grid on
xlabel('Velocity')
ylabel('Altitude')
zlabel('second')

set(gcf,'position',[x0,y0,width,height])
% title('a) Settling time','FontSize',12) 

% subplot(1,3,2)
figure
surf(SSE)
view(45,45)
hold on
% plot(mean(SSE),'-x','Color','r')
xlim([1,11]);
ylim([1,11]);
set(gca,'XTickLabel',Speeds,'YTickLabel',Altitude,'Ytick',1:11,'XTick',1:11)
colorbar("north")

% set(gca,'XTickLabel',Speeds)
grid on
xlabel('Velocity')
ylabel('Altitude')
zlabel('RMS')
% title('b) Steady state error','FontSize',12) 
set(gcf,'position',[x0,y0,width,height])

% subplot(1,3,3)
figure()
surf(Ctrl_eff)
view(45,45)
hold on
xlim([1,11]);
ylim([1,11]);
% plot(mean(Ctrl_eff),'-x','Color','r')
set(gca,'XTickLabel',Speeds,'YTickLabel',Altitude,'Ytick',1:11,'XTick',1:11)
colorbar("north")
grid on

xlabel('Velocity')
ylabel('Altitude')
zlabel('RMS')
% title('c) Control effort','FontSize',12) 
set(gcf,'position',[x0,y0,width,height])

%% Bode plot
figure
bode(SysLICg0(1,1,:,:),SysLICg(1,1,:,:),'--k');
legend('Sensivity function','Comprementary')
% b2.Responses(1).LineWidth = 1.25;
hold on 

%% inner loop
w = logspace(-4,4,1000);
for i=1:11
    for j=1:11
        [mag2,phase2,wout2] =bode(SysLICg0(2,1,i,j),w); % loop transfer
        subplot(2,1,1)
        semilogx(wout2, 20*log10(squeeze(mag2)),'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        xticklabels("")
        hold on

        subplot(2,1,2)
        semilogx(wout2, squeeze(phase2),'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        xlim([1e-3,1e3])
        xlabel('Frequency (rad/s)');
        hold on
    end
end

for i=1:11
    for j=1:11
        [mag3,phase3,wout3] =bode(1/(1+SysLICg0(2,1,i,j)),w); %sensivity
        subplot(2,1,1)
        semilogx(wout3, 20*log10(squeeze(mag3)),'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        xticklabels("")
        hold on

        subplot(2,1,2)
        semilogx(wout3, squeeze(phase3),'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        xlim([1e-3,1e3])
        xlabel('Frequency (rad/s)');
        hold on
    end
end

for i=1:11
    for j=1:11
        [mag,phase,wout] =bode(SysLICg(2,1,i,j),w); %transfer function
        % [mag2,phase2,wout2] =bode(SysLIC(1,1,i,j),w);
        s1=subplot(2,1,1);
        semilogx(wout, 20*log10(squeeze(mag)),'--k', 'LineWidth',1.25)
        hold on
        % semilogx(wout2, 20*log10(squeeze(mag2)),'--k', 'LineWidth',1.25)
        xlim([1e-3,1e3])
        ylabel('Magnitude (dB)');
        xticklabels("")
        
        grid
        s2=subplot(2,1,2);
        semilogx(wout, squeeze(phase),'--k', 'LineWidth',1.25)
        ylabel('Phase (deg)');
        hold on
        grid
    end
end
% legend('Sensivity function','Comprementary')


pos1 = get(s1, 'Position'); % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0];
set(s2, 'Position',new_pos1 );% set new position of current sub - plot


%%
figure() %% full loop
d2r=pi/180;
pid_h2gam=pid(0.2,0.005);lowpass=tf(1,[2,1]);

OL=d2r*pid_h2gam*lowpass*SysLICg(2,1,:,:);
OuterloopAllmargin=allmargin(OL);

Sen=1/(1+d2r*pid_h2gam*lowpass*SysLICg(2,1,:,:));


for i=1:11
    for j=1:11
        [mag5,phase5,wout5] =bode(OL(1,1,i,j),w); % loop transfer
        subplot(2,1,1)
        semilogx(wout5, 20*log10(squeeze(mag5)),'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        xticklabels("")
        hold on

        subplot(2,1,2)
        semilogx(wout5, squeeze(phase5),'Color',[0.3 0.3 0.3], 'LineWidth',1.25)
        xlim([1e-6,1e3])
        xlabel('Frequency (rad/s)');
        hold on
    end
end

for i=1:11
    for j=1:11
        [mag6,phase6,wout6] =bode(Sen(1,1,i,j),w); % Sensivity
        s1=subplot(2,1,1);
        semilogx(wout6, 20*log10(squeeze(mag6)),'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        xticklabels("")
        hold on

        s2=subplot(2,1,2);
        semilogx(wout6, squeeze(phase6),'Color',[0.5 0.5 0.5], 'LineWidth',1.25)
        xlim([1e-6,1e3])
        xlabel('Frequency (rad/s)');
        hold on
    end
end

w = logspace(-6,4,1000);
for i=1:11
    for j=1:11
        [mag4,phase4,wout4] =bode(CL(1,1,i,j),w); % transfer function
        subplot(2,1,1)
        semilogx(wout4, 20*log10(squeeze(mag4)),'--k', 'LineWidth',1.25)
        xticklabels("")
        hold on

        subplot(2,1,2)
        semilogx(wout4, squeeze(phase4),'--k', 'LineWidth',1.25)
        xlim([1e-6,1e3])
        xlabel('Frequency (rad/s)');
        hold on
    end
end



pos1 = get(s1, 'Position'); % gives the position of current sub-plot
new_pos1 = pos1 +[0 -0.37 0 0];
set(s2, 'Position',new_pos1 );% set new position of current sub - plot




