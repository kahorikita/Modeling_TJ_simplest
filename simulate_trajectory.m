clear
clc

delt = .001; % time step length in secs
plant.delt = delt;
Tmax = .25; % max. time of simulation in s 
imax = ceil(Tmax/plant.delt); % max timestep to compare between model and data
y = 1:1:1893;
Tmax_sim = (size(y,2))*plant.delt + .5; % always simulate 500 ms further than the data

% params for initial value of target jump
Hz = 130;
Tmax = .25; % max. time of simulation in s 
imax = ceil(Tmax/plant.delt); % max timestep to compare between model and data

Gstock = 0.14;
% Sstock = [0.001 0.002 0.003 0.004 0.005];
Sstock = [0.01 0.02 0.03 0.04 0.05];
Ustock = [0.1 0.2 0.3 0.4 0.5];

% change sigma
figure(1); clf; hold on;
for i = 1:5
    Ginit = Gstock;
    Sinit = Sstock(i);
    Uinit = Ustock(3);
    
    X = [Ginit Sinit Uinit];
    sim = sim_vel_X1D_BSL(X,plant,Tmax_sim);
    
    plot(delt*(1:imax),sim.acc(1:imax))
    legend({'sigma=0.01','sigma=0.02','sigma=0.03','sigma=0.04','sigma=0.05'},'Location','best');
end
ylim([0 70])
% xlim([0.1 0.2])

% change A
figure(2); clf; hold on;
for j = 1:5
    Ginit = Gstock;
    Sinit = Sstock(3);
    Uinit = Ustock(j);
    
    X = [Ginit Sinit Uinit];
    sim = sim_vel_X1D_BSL(X,plant,Tmax_sim);
    
    plot(delt*(1:imax),sim.acc(1:imax))
    legend({'A=0.1','A=0.2','A=0.3','A=0.4','A=0.5'},'Location','best');
end
ylim([0 70])
% xlim([0.12 0.16])

