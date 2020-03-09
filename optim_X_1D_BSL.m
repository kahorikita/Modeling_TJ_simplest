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
imax = ceil(Tmax/plant.delt); % max timestep to compare between model and data

% params for initial value of target jump
Hz = 130;

% G = timing of target jump, U = amplitude 
% Gstock = [0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15];
% Ustock = [0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.011];

% Gstock = 0.05;  % 0.07 mean vel for all
% Sstock = 0.002; % 0.002 mean vel for all
% Ustock = 0.005; % 0.005 mean vel for all
Gstock = [0.05 0.06 0.07 0.08 0.09];
Sstock = [0.001 0.002 0.003 0.004 0.005];
Ustock = [0.004 0.005 0.006 0.007 0.008];

% c = 5;
% fhandle = figure(c); clf; hold on
% set(fhandle, 'Position', [200, 100, 900, 650]); % set size and loction on screen
% set(fhandle, 'Color','w') % set background color to white 
% set(gca,'FontSize',10);

c = 1;

Xopt = [];


% switch data used in simulation
% simulate velocity / mean of all subs / bootstrap for each sub
disp('1:simulated vel, 2:mean of all subs, 3:bootstrap acc data, 4:mean of each sub, 5:mean acc of all subs, 6:mean acc of each sub');  
m = input('Choose data: ');
switch m
    case 1 % to simulate velocity profiles for parameter recovery
        disp('simulated vel')  
%         TrueParam = Result{1,1};
        TrueParams = [.15 .025 1.56 0.6 0.06 -0.87]; %[TJ, 2 params of gausian distribution, 4 params of L] 
        %for numofsim = 1:100     
            sim = sim_vel_X1D_BSL(TrueParams,H,plant,Tmax);
        %end
        y = sim.convo;
        numofbootstraps = 1;
        numofsubjects = 1;
    case 2 % average velocity of all subjects (experiment data)
        disp('mean vel of all subs')
        load allsub_meanvel.mat
        y_130 = mean_vel{1}(3,:)-mean_vel{1}(4,:); % target on X, jump +/-y, BSL
        y = resample(y_130,1/delt,Hz)*Hz; % resample to 1000Hz
        y = y(101:end); % subtract 100 ms instrument delay
        y = y - mean(y(1:100)); % subtract baseline
        numofbootstraps = 1;
        numofsubjects = 1;
    case 3 % each subject acc data with bootstrap
        disp('bootstrap acc data')
        load AllData.mat
        numofsubjects = 20;
        numofbootstraps = 5;
    case 4 % mean velocity of each subject
        disp('mean vel of each sub')
        load AllData.mat
        numofsubjects = 20;
        numofbootstraps = 1;
    case 5 % mean acc of all subjects
        disp('mean acc of all subs')
        load allsub_meanvel.mat        
        y_130 = mean_vel{1}(3,:)-mean_vel{1}(4,:); % target on X, jump +/-y, BSL
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
        numofsubjects = 20;
        numofbootstraps = 1;
        bs = 1;
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


for subject = 1:numofsubjects
% for subject = 1:1

    if m == 3
        datap = data(subject,3+2).Vel_CrX_post*Hz; % target jump +y data{1,3} data(subject,1:2)=TR
        datam = data(subject,4+2).Vel_CrX_post*Hz; % target jump -y data{1,4}
        ybs = get_bootstrap(datap,datam,numofbootstraps,delt,Hz);
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
    
    
    for bs = 1:numofbootstraps % number of bootstrap     
        count = 1;
    
        if m == 3
            % set one bootstrap data, acc
            y = diff(ybs(bs,:))/delt;
        end
        
        for i = 1:length(Gstock)    % istart:istart+numofsteps
                
            Ginit = Gstock(i);
            
            for j = 1:length(Sstock)
                
                Sinit = Sstock(j);
                
                for k = 1:length(Ustock)
                    
                Uinit = Ustock(k);
                Xinit = [Ginit Sinit Uinit];
                f_targ = @(X) get_error_X1D_BSL(X,y,plant,Tmax,m);
                
                switch optm
                    case 1 % fmincon
                        Aeq = [];
                        beq = [];
                        % used for TA and myu
%                         lb =  [0.05   0]; % Lower bounds
%                         ub =  [0.5 Inf]; % Upper bounds % 0.3 
                        lb = [];
                        ub = [];
                        
                        Xopt = fmincon(f_targ,Xinit,[],[],Aeq,beq,lb,ub);   
                        
                    case 2 % bads
                        lb =  [0.04 0.0001 0.001]; % Lower bounds
                        plb = [0.07 0.0001 0.001]; % Plausible Lower bounds
                        pub = [0.2 0.1 0.1]; % Plausible Upper bounds
                        ub =  [0.3 0.2 0.3]; % Upper bounds % 0.3
                        
                        if m == 3 || m == 6 % 3:bootstrap acc, 6:mean acc of each sub
                            temp = bads(f_targ,Xinit,lb,ub,plb,pub);
                            Tmax_sim = (size(y,2))*plant.delt + .5;
                            opt = sim_vel_X1D_BSL(temp,plant,Tmax_sim);
                            error = nanmean((y(1:imax-1)-opt.acc(1:imax-1)).^2); % error
                            Xopt{subject,bs}(count,:) = [Xinit temp error];     
                        else
                            Xopt = bads(f_targ,Xinit,lb,ub,plb,pub);    
                        end
                            
                    otherwise
                end
                
                count = count + 1;
                
                end
            end
        end
      
        Tmax_sim = (size(y,2))*plant.delt + .5;
        if m == 3 % bootstrap, acc
            [M(subject,bs),I(subject,bs)] = min(Xopt{subject,bs}(:,7)); % choose G and U with minimum error
            opt = sim_vel_X1D_BSL(Xopt{subject,bs}(I(subject,bs),4:6),plant,Tmax_sim); 
            
            figure(10+subject); 
            subplot(1,5,bs); hold on;
            plot(delt*(1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            plot(delt*(1:imax),opt.acc(:,1:imax),'color',cols(1,:,c),'linewidth',1.5)
            plot([Xopt{subject,bs}(I(subject,bs),4),Xopt{subject,bs}(I(subject,bs),4)],[-300,400],'k','linewidth',1.5);
%             legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Acceleration (m/s*s)','FontSize',10)
            xlim([0 1]);
                   
        elseif m == 6 % input data = accerelation, each sub
            [M(subject),I(subject)] = min(Xopt{subject}(:,7)); % choose G and U with minimum error 
            opt = sim_vel_X1D_BSL(Xopt{subject}(I(subject),4:6),plant,Tmax_sim); 

            figure(6); hold on;
            subplot(4,5,subject); hold on;
            plot(delt*(1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            plot(delt*(1:imax),opt.acc(:,1:imax),'color',cols(1,:,c),'linewidth',1.5)
            plot([Xopt{subject}(I(subject),4),Xopt{subject}(I(subject),4)],[-300,400],'k','linewidth',1.5);
%             legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Acceleration (m/s*s)','FontSize',10)
            xlim([0 1]);
            
        elseif m == 5 % input data = accerelation, mean of all subs
            Tmax_sim = (size(y,2))*plant.delt + .5;
            opt = sim_vel_X1D_BSL(Xopt,plant,Tmax_sim);    
            
            figure(7);
            subplot(1,2,1); hold on;
            plot(delt*(1:length(dvel)),dvel,'color',cols(4,:,c),'linewidth',1.5)
            plot(delt*(1:imax),opt.convo(:,1:imax),'color',cols(1,:,c),'linewidth',1.5)
            plot([Xopt(1),Xopt(1)],[-0.2,0.6],'k');
            legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Velocty (m/s)','FontSize',10)
            
            subplot(1,2,2); hold on;
            plot(delt*(1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            plot(delt*(1:imax),opt.acc(:,1:imax),'color',cols(1,:,c),'linewidth',1.5)
            plot([Xopt(1),Xopt(1)],[-300,400],'k');
            legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Acceleration (m/s*s)','FontSize',10)
        else
            figure(1); clf; hold on;
            opt = sim_vel_X1D_BSL(Xopt,plant,Tmax_sim);    
        
            plot(delt*(1:length(y)),y,'color',cols(4,:,c),'linewidth',1.5)
            plot(delt*(1:imax),opt.x(:,1:imax),'color',cols(1,:,c),'linewidth',1.5)
%             legend('simulated data','fit')
            xlabel('Time (s)','FontSize',10)
            ylabel('Velocty (m/s)','FontSize',10)
        end
        
        
    end
end

% save 'TJModelFits.mat' Resultopt Result_all


toc
