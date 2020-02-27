function sim = sim_vel_X1D_BSL(X,plant,Tmax)

delt = plant.delt;

Tjump = ceil(X(1)/delt); % mean time of Target Jump is gittered (gausian distribition)
uamp = X(2); % amplitude of fitting function

% adjust high level parameters
start = 0; % starting position of the hand
nstep = ceil(Tmax/delt); % number of timesteps to simulate

% initialize the state vector, x
x = 1:nstep;
x(1,1) = start; % set starting hand position

% simulate velocity with parameters X(1), X(2) 
Y = uamp*power(max(x-Tjump,0),2)*delt;

% save data to sim
sim.x = Y;
sim.acc = diff(Y)/delt;
sim.T = Tjump;
sim.delt = delt;
sim.plant = plant;
sim.X = X;

% figure(tr); cla
% hold on;
% % xax1 = 1:length(sim.x(2,:));
% xax1 = 1:length(sim.convo);
% time1 = 7.7*xax1;
% % plot(time1,sim.x(2,:),'g')
% plot(time1,sim.convo,'g')
% xax2 = 1:length(y);
% time2 = 7.7*xax2;
% plot(time2,y,'r')
% 
% plot([7.7*T(1) 7.7*T(1)],[-0.1 0.4],'k');
% plot([7.7*T(2) 7.7*T(2)],[-0.1 0.4],'b');

% axis([0 1400 -0.05 0.25]);




