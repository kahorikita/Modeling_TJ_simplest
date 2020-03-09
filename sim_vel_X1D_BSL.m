function sim = sim_vel_X1D_BSL(X,plant,Tmax)

delt = plant.delt;

Tjump = ceil(X(1)/delt); % mean time of Target Jump is gittered (gausian distribition)
sigma = abs(X(2))/delt; % variance of response time 
uamp = X(3); % amplitude of fitting function
% uamp = X(2); % amplitude of fitting function

% adjust high level parameters
start = 0; % starting position of the hand
nstep = ceil(Tmax/delt); % number of timesteps to simulate

% initialize the state vector, x
x = 1:nstep;
x(1,1) = start; % set starting hand position

% simulate velocity with parameters X(1), X(2) 
Y = uamp*power(max(x-Tjump,0),2)*delt;

% convolve with Gaussian to account for variable response onset across
% trials
t_min = ceil(sigma)*3; % set convolution range to 3 standard deviations
y_norm = normpdf(-t_min:t_min,0,sigma); % gausian distribution

% save data to sim
sim.x = Y;
sim.convo = conv(Y,y_norm,'same');
% sim.acc = diff(Y)/delt;
sim.acc = diff(sim.convo)/delt;
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




