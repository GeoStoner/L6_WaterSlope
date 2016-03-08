% Ryan Stoner. March 7, 2016 for modeling in the Earth Sciences
clear
clc
%% initialize
% Rainwater movement
R = 0.005;              % m/s, recharge rate
I = 0.0045;              % m/s, infiltration rate

% Creating the basement
zmax = 8;              % m, initial height
s = 0.05;               % slope
xmax = 100;             % m
xmin = 0;               % m
dx = 2;                 % m
x = xmin:dx:xmax;       % m
zbas = zmax - s*x;         % m

% initial water height
H = zeros(1,length(x)); % m

% initializing time
tmax = 30000;             % s
dt = 0.1;                 % s
t = 0:dt:tmax;          % s
nplots = 20;
tplot = tmax/nplots;

nsteps = length(t);
nframe = 0;

% initializing constants

n = 0.030;              % roughness coefficient
se = s;                 % energy slope

%% Loop

for i=1:nsteps

 ubar = (1/n)*H.^(2/3)*se^(1/2);
 q = ubar.* H;
 dqdx = diff(q)/dx;
 dqdx = [dqdx 0];
 dhdt = -dqdx + R - I;

 % Update water height
 H = H + dhdt* dt;
 z = zbas + H;
 
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

 end

end

%% Close

% Analytical Solution

Han = nthroot((((R-I)*x)/((1/n)*se^1/2)).^3,5);
