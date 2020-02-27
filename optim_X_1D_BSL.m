%%% Linit is calcuted by optim_L_1D

%%% -define the true parameters
%%% -simulate model to get 'true' trajectory
%%% -save trajectory to compare with candidate model fits in get_error function
clear
clc

tic

cols(:,:,1) = [ 0 210 255; 255 210 0; 0 0 0; 210 0 255]/256;
cols(:,:,2) = [ 0 155 255; 255 100 0; 0 0 0; 155 0 255]/256;
cols(:,:,3) = [ 0 100 255; 255 0 0; 0 0 0; 100 0 255]/256;

% for simulate velocity 
delt = .001; % time step length in secs
plant.delt = delt;

Tmax = .25; % max. time of simulation in s
len = ceil(Tmax/delt);

% params for initial value of target jump
% istart = 30; %20
% numofsteps = 1; %15
Hz = 130;


%%

% switch data used in simulation
% simulate velocity / mean of all subs / bootstrap for each sub
disp('1:simulated vel, 2:mean of all subs, 3:bootstrap data, 4:mean of each sub, 5:mean acc of all subs, 6:mean acc of each sub');  
m = input('Choose data: ');
switch m
    case 1 % to simulate velocity profiles for parameter recovery
        disp('simulated vel')  
%         TrueParam = Result{1,1};
        TrueParams = [.15 .025 1.56 0.6 0.06 -0.87]; %[TJ, 2 params of gausian distribution, 4 params of L] 
        %for numofsim = 1:100 % it's different from the "tr" below       
            sim = sim_vel_X1D_BSL(TrueParams,H,plant,Tmax)
        %end
        y = sim.convo;
        numofbootstraps = 1;
        numofsubjects = 1;
        
    case 2 % average velocity of all subjects (experiment data)
        disp('mean vel of all subs')
        load allsub_meanvel.mat
%         y = mean_vel{2}(3,1:len)*Hz-mean_vel{2}(4,1:len)*Hz; % target on X, jump +/-y, MR
        y_130 = mean_vel{1}(3,:)-mean_vel{1}(4,:); % target on X, jump +/-y, BSL
        y = resample(y_130,1/delt,Hz)*Hz; % resample to 1000Hz
        y = y(101:end); % subtract 100 ms instrument delay
        y = y - mean(y(1:100)); % subtract baseline
        numofbootstraps = 1;
        numofsubjects = 1;
    case 3 % each subject data with bootstrap
        disp('bootstrap data')
        load AllData.mat
        numofsubjects = 20;
    case 4 % mean velocity of each subject
        disp('mean vel of each sub')
        load AllData.mat
        numofsubjects = 20;
        numofbootstraps = 1;
    case 5 % mean acc of all subjects
        disp('mean acc of all subs')
        load allsub_meanvel.mat        
        y_130 = mean_vel{1}(3,:)-mean_vel{1}(4,:); % target on X, jump +/-y, BSL
%         y = mean_vel{2}(3,1:len)/(1/Hz)-mean_vel{2}(4,1:len)/(1/Hz); % target on X, jump +/-y, MR
        y = resample(y_130,1/delt,Hz)*Hz; % resample to 1000Hz
        y = y(101:end); % subtract 100 ms instrument delay
        y = y - mean(y(1:100)); % subtract baseline
        dvel = y;
        y = diff(y)/delt;
        numofbootstraps = 1;
        numofsubjects = 1;
    case 6 % mean acc of each subject
        disp('mean acc of each sub')
        load AllData.mat
        numofsubjects = 2;
        numofbootstraps = 1;  
    otherwise
        disp('other value')
end

disp('1:MLE, 2:BADS');
optm = input('Choose an optimization method: '); %
switch optm
    case 1
        disp('MLE')
    case 2
        disp('BADS')
end

c = 5;
fhandle = figure(c); clf; hold on
set(fhandle, 'Position', [200, 100, 900, 650]); % set size and loction on screen
set(fhandle, 'Color','w') % set background color to white 
set(gca,'FontSize',10);

c = 6;
fhandle = figure(c); clf; hold on
set(fhandle, 'Position', [300, 100, 900, 650]); % set size and loction on screen
set(fhandle, 'Color','w') % set background color to white 
set(gca,'FontSize',10);

c = 1;

% Gstock = [0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15];
% Ustock = [0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.011];
Gstock = [0.16 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15];
Ustock = [0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05];

Xopt = [];

for subject = 1:numofsubjects

    if m == 3
        datap = data(subject,3+2).Vel_CrX_post*Hz; % target jump +y data{1,3} data(subject,1:2)=TR
        datam = data(subject,4+2).Vel_CrX_post*Hz; % target jump -y data{1,4}
        ybs = get_bootstrap(datap,datam,len,numofbootstraps);
    elseif m == 4 % mean vel of each subject         
        y_130 = nanmean(data(subject,3+2).Vel_CrX_post)-nanmean(data(subject,4+2).Vel_CrX_post);
        y = resample(y_130,1/delt,Hz)*Hz; % resample to 1000Hz
        y = y(101:end); % subtract 100 ms instrument delay
        y = y - mean(y(1:100)); % subtract baseline
    elseif m == 6 % mean acc of each subject
        y_130 = nanmean(data(subject,3+2).Vel_CrX_post)-nanmean(data(subject,4+2).Vel_CrX_post);
        y = resample(y_130,1/delt,Hz)*Hz; % resample to 1000Hz
        y = y(101:end); % subtract 100 ms instrument delay
        y = y - mean(y(1:100)); % subtract baseline
        dvel = y;
        y = diff(y)/delt;
    end
    
    for tr = 1:numofbootstraps % number of bootstrap     
        
        if m == 3
            % set one bootstrap data
            y = ybs(tr,:);
        end
        
        for i = 1:length(Gstock)    % istart:istart+numofsteps
                
            Ginit = Gstock(i);
            count = length(Gstock)*(i-1);
            
            for j = 1:length(Ustock)
                
                Uinit = Ustock(j);
                
%                 Ginit = .1; % initial guess for Gaussian parameters % ideal for mean all subs 0.1
%                 Uinit = .006; % initial guess for u amplitude % ideal for mean all subs 0.006
                
                Xinit = [Ginit Uinit];
                
                f_targ = @(X) get_error_X1D_BSL(X,y,plant,Tmax,m);
                
                switch optm
                    case 1 % fmincon
                        Aeq = [];
                        beq = [];
                        % used for TA and myu
                        lb = [0 0.001 -2 -2 -2 -2]; % Lower bounds
                        ub = [0.5 0.5 2 2 2 2]; % Upper bounds
                        
                        Xopt = fmincon(f_targ,Xinit,[],[],Aeq,beq,lb,ub);
                        
                    case 2 % bads
                        lb =  [0.05   0]; % Lower bounds
                        plb = [0.07 0]; % Plausible Lower bounds
                        pub = [0.2 1]; % Plausible Upper bounds
                        ub =  [0.3 Inf]; % Upper bounds % 0.3       
                        
                        if m == 6
                            Xopt{subject}(count+j,1) = Ginit;
                            Xopt{subject}(count+j,2) = Uinit;
                            Xopt{subject}(count+j,3:4) = bads(f_targ,Xinit,lb,ub,plb,pub);
                            opt = sim_vel_X1D_BSL(Xopt{subject}(count+j,3:4),plant,Tmax);  
                            Xopt{subject}(count+j,5) = nanmean((y(1:len)-opt.x(1:len)).^2);
                        else
                            Xopt = bads(f_targ,Xinit,lb,ub,plb,pub);    
                        end
                            
                    otherwise
                end
                
            end
        end
      
        if m == 6 % input data = accerelation, each sub
            figure(15+subject); clf; hold on; 
            time = 0.001;
%             Xopt_org = Xopt;
%             Iz = find(Xopt{subject}(:,1)==0);
%             Xopt{subject}(Iz,:) = [];
            [M(subject),I(subject)] = min(Xopt{subject}(:,5)); % choose G and U with minimum error 
            opt = sim_vel_X1D_BSL(Xopt{subject}(I(subject),3:4),plant,Tmax); 
%             subplot(4,5,subject); hold on;
            plot(time*(1:length(dvel)),dvel,'color',cols(4,:,c),'linewidth',1.5)
            plot(time*(1:len),opt.x(:,1:len),'color',cols(1,:,c),'linewidth',1.5)
%             plot(time*(1:len),cumsum(opt.acc(:,1:len))/1000,'r','linewidth',1.5)
            plot([Xopt{subject}(I(subject),3),Xopt{subject}(I(subject),3)],[-0.1,0.3],'k','linewidth',1.5);
            legend('data','model','response time')
            xlabel('Time (s)','FontSize',10)
            ylabel('Velocty (m/s)','FontSize',10)
            xlim([0 1]);
            
            figure(6); hold on;
            subplot(4,5,subject); hold on;
            plot(time*(1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            plot(time*(1:len-1),opt.acc(:,1:len-1),'color',cols(1,:,c),'linewidth',1.5)
            plot([Xopt{subject}(I(subject),3),Xopt{subject}(I(subject),3)],[-1,3],'k','linewidth',1.5);
%             legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Acceleration (m/s*s)','FontSize',10)
            xlim([0 1]);
            
        elseif m == 5 % input data = accerelation, mean of all subs
            opt = sim_vel_X1D_BSL(Xopt,plant,Tmax);    
            figure(subject+10); clf; hold on;
            time = 0.001; 
            subplot(1,2,1); hold on;
            plot(time*(1:length(dvel)),dvel,'color',cols(4,:,c),'linewidth',1.5)
            plot(time*(1:len),opt.x(:,1:len),'color',cols(1,:,c),'linewidth',1.5)
%             plot(time*(1:len),cumsum(opt.acc(:,1:len))/1000,'r','linewidth',1.5)
            plot([Xopt(1),Xopt(1)],[-0.2,0.6],'k');
            legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Velocty (m/s)','FontSize',10)
            
            subplot(1,2,2); hold on;
            plot(time*(1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            plot(time*(1:len-1),opt.acc(:,1:len-1),'color',cols(1,:,c),'linewidth',1.5)
            plot([Xopt(1),Xopt(1)],[-3,4],'k');
            legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Acceleration (m/s*s)','FontSize',10)
        else
            opt = sim_vel_X1D_BSL(Xopt,plant,Tmax);    
            figure(1); clf; hold on;
            time = 0.001;
            figure(1); clf; hold on;
            plot(time*(1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            plot(time*(1:len),opt.x(:,1:len),'color',cols(1,:,c),'linewidth',1.5)
%             plot(population2);
%             plot((1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Velocty (m/s)','FontSize',10)
        end
        
        
    end
end

% save 'TJModelFits.mat' Resultopt Result_all


toc
