function y = get_bootstrap(datap,datam,numofbootstrap,delt,Hz)

for i=1:numofbootstrap
    % 1. bootstrap and calculate mean velocity
    % for target jump +y
    % datap = data{1,3}.Vel_CrX_post*Hz;
    bootp = datasample(datap,size(datap,1));
    
    % 2. for target jump -y
    % datam = data{1,4}.Vel_CrX_post*Hz;
    bootm = datasample(datam,size(datap,1));
    
    % 3. calculate diff_vel of mean velocities
    ytmp_130 = nanmean(bootp) - nanmean(bootm);
    ytmp = resample(ytmp_130,1/delt,Hz)*Hz; % resample to 1000Hz
    ytmp = ytmp(101:end); % subtract 100 ms instrument delay
    ytmp = ytmp - mean(ytmp(1:100)); % subtract baseline
    y(i,:) = ytmp;
    
end

% 4. optimization

% 5. repeat 1-4

