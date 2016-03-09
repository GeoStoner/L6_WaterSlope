% Ryan Stoner. March 7, 2016 for modeling in the Earth Sciences
clear
clc
%% initialize
% Rainwater movement
R = 0.05;                         % m/s, recharge rate
I = 0.02;                         % m/s, infiltration rate

% Creating the basement
zmax = 8;                           % m, initial height
s = 0.05;                           % slope
xmax = 100;                         % m
xmin = 0;                           % m
dx = 1;                             % m
x = xmin:dx:xmax;                   % m
N = length(x);                      % used for matrix sizes
zbas = zmax - s*x;                  % m

% initial water height
initH =ones(1,N);                   % m
H = 0*initH;                        % m

% find edge values
hedge(1:N-1) = H(1:N-1)+diff(H)/2;

% initializing time
tmax = 25;                         % s
dt = 0.002;                        % s
t = 0:dt:tmax;                      % s
nplots = 60;                        % number of plots
tplot = tmax/nplots;                % time interval between plots

nsteps = length(t);                 % number of steps in loop
nframe = 0;                         % initializing counter for plotting

% initializing constants and looped values

n = 0.030;                          % roughness coefficient, gravel bed
se = s;                             % energy slope

ubar = zeros(1,N);
Q = zeros(1,N);
dQdx = zeros(1,N);
dhdt = zeros(1,N);

%% Loop

for i=1:nsteps
 
 % find mean speed of water/fluid
 ubar = (1/n)*hedge.^(2/3)*se^(1/2);
 
 % top of slope, so no water added from above, boundary condition
 Q(1) = 0; 
 Q(2:N) = ubar.* hedge;
 
 % find change in flux, add boundary condition to let water out of system
 dQdx(1:N-1) = diff(Q)/dx;
 dQdx(N) = dQdx(N-1);
 
 dhdt = -dQdx + R - I;

 % Update water height and hedge
 H = H + dhdt* dt;
 z = zbas + H;
 hedge(1:N-1) = H(1:N-1)+diff(H)/2;
 
 Hbelow = find(H<=0); % just like with corals, make sure water not negative
 H(Hbelow)=0;
 % Plotting a limited number of plots

 if(mod(t(i),tplot)==0)
 nframe = nframe+1;
    figure(1)
   
    plot(x,zbas,'k')
    
    % Filling in area
    baseval1 = -10000;
    
    
    % Plot of water height
    hold on
    plot(x,z,'r')
    
    baseval1 = -10000;
    a2 = area(x,z);
    ec = a2.FaceColor;
    a2.FaceColor = [0 0 1];
    baseval1 = -10000;
    a1 = area(zbas);
    ec = a1.FaceColor;
    a1.FaceColor = [0 1 1];
    % Converting time to a string to print string out

    trounded = round(t(i),0);
    strtime = num2str(trounded);


    txt = uicontrol('Style','text',...
      'Position',[430 0 120 20],...
      'Fontsize',12,...
      'FontWeight','Bold',... 
      'HorizontalAlignment','left',...
      'String',['Time: ' strtime ' seconds']);

    xlabel('distance (km)')
    ylabel('water height (m)')
    xlim([xmin+dx xmax])
    ylim([zmax-s*x(N) zmax])
    title('Movement of Water over Land Surface')


    pause(0.1)
    hold off
 end

end

%% Close

% Analytical Solution
Qan = (R-I)*x;
Han = nthroot(( ((R-I)*x*n)/(se^(1/2))).^3,5);

% Printing verification of numerical solution and analytical solution
 if(abs(Qan-Q)<=0.001)
     outstrq = 'The analytical and numerical solution for Q match \n'; 
     fprintf(outstrq);
 else
     outstrqneg = 'The analytical and numerical solution for Q DO NOT match \n'; 
     fprintf(outstrqneg);
 end
 
 % H is harder to match exactly, but smaller space and time steps improve
 % results.
  if(abs(Han-H)<=0.1)
     outstrH = 'The analytical and numerical solution for H match \n'; 
     fprintf(outstrH);
 else
     outstrHneg = 'The analytical and numerical solution for H DO NOT match \n'; 
     fprintf(outstrHneg);
 end
 

