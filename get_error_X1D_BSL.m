function [e2,A] = get_error_X1D_BSL(X,y,plant,Tmax,m)

Tmax_sim = (size(y,2))*plant.delt + .5; % always simulate 500 ms further than the data

% perform simulations
sim = sim_vel_X1D_BSL(X,plant,Tmax_sim); 

imax = ceil(Tmax/plant.delt); % max timestep to compare between model and data
if(imax>size(y,2))
    error('Tmax too long')
end

length_sim = size(sim.x,2);
length_data = size(y,2);

time = 0.001;

% calculate error between simulation and acual data 
if m == 5 || m == 6 || m == 3 % input is acceleration
    
    A = y(1:imax)/sim.acc(1:imax);
    e2 = nanmean((y(1:imax)-A*sim.acc(1:imax)).^2);
%     figure(3); clf; hold on
%     
%     plot(time*(1:length(y(1:imax))),y(1:imax),'m')
%     plot(time*(1:length(sim.acc(1:imax))),sim.acc(1:imax),'b')
%     plot(time*(1:length((y(1:imax)-sim.acc(1:imax)))),(y(1:imax)-sim.acc(1:imax)),'r')

%     plot(y,'m')
%     plot(sim.acc,'b')
%     plot((y(1:imax)-sim.acc(1:imax)),'r')

else % input is velocity
%     e2 = nanmean((y(1:imax)-sim.x(1:imax)).^2);
    A = y(1:imax)/sim.convo(1:imax);
    e2 = nanmean((y(1:imax)-A*sim.convo(1:imax)).^2);
    figure(3); clf; hold on

    plot(time*(1:length(y)),y,'m')
    plot(time*(1:length(A*sim.convo)),sim.convo,'b')
    plot(time*(1:length((y(1:imax)-A*sim.convo(1:imax)))),(y(1:imax)-A*sim.convo(1:imax)),'r')

%     plot(y,'m')
%     plot(sim.x,'b')
%     plot((y(1:imax)-sim.x(1:imax)),'r')
end

sim.e2 = e2;

% disp(num2str(e2))
%
% figure(3); clf; hold on
% plot(y,'m')
% plot(sim.convo,'b')
% plot((y(1:imax)-sim.convo(1:imax)),'r')
%}

