%% Example ASA
% Find minima of a function with Simulated Annealing
% Manel Soria, David de la Torre and Arnau Miro - ETSEIAT

% Clean-up
close all;
clear;

% Our test is a R^2->R function based on Rastrigin function. It is 
% challenging because it has infinite local extrema, located at integer
% numbers (ie, 8,-9) 
% The global minimum is at (1,1), and its value is 0
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define SA parameters
nitemax = 100; % Maximum number of iterations
mu = 0.2; % Thermal transition probability parameter
goal = 1E-5; % Target fitness value

% Define SA info data
ninfo = 10; % Verbosity level (print every # iterations)
label = 10; % Label (identification purposes)

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define SA functions
fitfun = @(x) ras(x(1),x(2)); % Fitness function - TO BE MINIMIZED
mutfun = @(x,f) x + ranrange(-0.1,0.1,2); % Mutation: small random mov
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Initial guess
A0 = [0.23; 1.45];

% Execute Simulated Annealing
[bestIndSA,bestFitSA,...
    nite,history] = asa(label,ninfo,...
                        A0,nitemax,mu,goal,...
                        fitfun,mutfun,prifun);

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1e-8,'Display','none');
[bestIndFMS,bestFitFMS] = fminsearch(fitfun,bestIndSA,options);

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('ASA \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndSA,bestFitSA);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndFMS,bestFitFMS);

%% Plot history

% Create figure
fh = figure();
fh.Position = [400,200,900,600];

% Plot history
semilogy(history,'o-');

% Beautify plot
grid minor;
title('Rastrigin function | Simulated Annealing optimization');
xlabel('Generation [#]');
ylabel('Fitness function value');

