% Ryan Stoner. March 7, 2016 for modeling in the Earth Sciences
clear
clc
%% initialize
% Rainwater movement
R = 0.005;                          % m/s, recharge rate
I = 0.0045;                         % m/s, infiltration rate

% Creating the basement
zmax = 8;                           % m, initial height
s = 0.05;                           % slope
xmax = 100;                         % m
xmin = 0;                           % m
dx = 2;                             % m
x = xmin:dx:xmax;                   % m
N = length(x);                      % used for matrix sizes
zbas = zmax - s*x;                  % m

% initial water height
initH =ones(1,N);           % m
H = 0*initH;                        % m

% find edge values
hedge(1:N-1) = H(1:N)-1)+diff(H)/2;

% initializing time
tmax = 30000;                       % s
dt = 0.1;                           % s
t = 0:dt:tmax;                      % s
nplots = 20;                        % number of plots
tplot = tmax/nplots;                % time interval between plots

nsteps = length(t);                 % number of steps in loop
nframe = 0;                         % initializing counter for plotting

% initializing constants

n = 0.030;                          % roughness coefficient, gravel bed
se = s;                             % energy slope

%% Loop

for i=1:nsteps
 
 % find mean speed of water/fluid
 ubar = (1/n)*hedge.^(2/3)*se^(1/2);
 q(2:N) = ubar.* hedge;
 
 % top of slope, so no water added from above, boundary condition
 q(1) = 0; 
 
 % find change in flux, add boundary condition to let water out of system
 dqdx(1:N-1) = diff(q)/dx;
 dqdx(N) = dqdx(N-1);
 
 dhdt = -dqdx + R - I;

 % Update water height
 H = H + dhdt* dt;
 z = zbas + H;
 
 Hbelow = find(H<0); % just like with corals, we make sure water not negative
 H(Hbelow)=0;
% Plotting a limited number of plots

 if(mod(t(i),tplot)==0)
nframe = nframe+1;
    figure(1)
   
    plot(x,zbas,'k')
    hold on
    plot(x,z,'r')

% Converting time to a string to print string out

    trounded = round(t(i),0);
    strtime = num2str(trounded);


    txt = uicontrol('Style','text',...
      'Position',[430 0 80 20],...
      'Fontsize',12,...
      'FontWeight','Bold',... 
      'HorizontalAlignment','left',...
      'String',['Time: ' strtime]);

    xlabel('distance (km)')
    ylabel('elevation (m)')
    title('Movement of Water over Land Surface')


    pause(0.1)
    hold off
 end

end

%% Close

% Analytical Solution

Han = nthroot((((R-I)*x)/((1/n)*se^1/2)).^3,5);
