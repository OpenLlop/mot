%% Example AGA
% Find minima of a function with a Genetic Algorithm
% Manel Soria, David de la Torre and Arnau Miro - ETSEIAT

% Clean-up
close all;
clear;

% Our test is a R^2->R function based on Rastrigin function. It is 
% challenging because it has infinite local extrema, located at integer
% numbers (ie, 8,-9) 
% The global minimum is at (1,1), and its value is 0
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define GA parameters
ng = 50; % Number of generations
np = 200; % Population size
N = [3,... % Number of elites
    floor(np*0.1),... % Number of mutants
    floor(np*0.05),...% Number of newcommers
    floor(np*0.2)]; % Number of parents
goal = 1E-5; % Target fitness value

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define GA functions
funique = @(x) x; % Discard identical individuals: currently not in use
fitfun = @(x) ras(x(1),x(2));  % Fitness function - TO BE MINIMIZED
mutfun = @(x,f) x + ranrange(-0.1,0.1,2); % Mutation: small random mov
repfun = @(x,y) (x+y)/2; % Reproduction: average
ranfun = @() ranrange(-200,200,2); % Random individual
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repetability in the GA

% Define GA info data
info = 2; % Verbosity level
label = 10; % Label

% Execute Genetic Algorithm
[bestPopGA,bestFitValGA,nite,history] = aga(info,label, ...
                 np,...
                 ng,N,goal,...
                 funique,fitfun,mutfun,repfun,ranfun,prifun);
    % pop2 is the sorted population after iterating; pop2{1} is the best solution available
    % best = fitfun(pop2{1}) is the optimal value
    % nite is the number of iterations performed
bestIndGA = bestPopGA{1}; % Best Individual of GA algorithm

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1e-8,'Display','none');
[bestIndFMS,bestFitValFMS] = fminsearch(fitfun,bestIndGA,options);

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('AGA \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndGA,bestFitValGA);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndFMS,bestFitValFMS);

%% Plot history

% Create figure
fh = figure();
fh.Position = [400,200,900,600];

% Plot history
semilogy(history,'o-');

% Beautify plot
grid minor;
title('Rastrigin function | Genetic Algorithm optimization');
xlabel('Generation [#]');
ylabel('Fitness function value');

