%% Example AGA
% Find minima of a function with Genetic Algorithm + Pareto fronts (GAP)
%
% Programmers:   Manel Soria         (UPC/ETSEIAT)
%                David de la Torre   (UPC/ETSEIAT)
% Date:          23/11/2023
% Revision:      1
close all; clc;

%% AGA

% Define heuristic function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 1; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})

% Define AGA algorithm parameters
goal = -1; % Target fitness value
npf = 20; % Maximum number of pareto fronts
ng = 10; % Number of generations
np = 500; % Population size
N = [1,... % Number of elites
    floor(np*0.1),... % Number of mutants
    floor(np*0.05),...% Number of newcomers
    floor(np*0.2)]; % Number of parents

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a and b
perfun = @(x) 2 * (x(1) + x(2)); % Perimeter

% Equilateral triangle
x0 = [5,5]'; % Offset
xa = x0 + [0,0]';
xb = x0 + [2,0]';
xc = x0 + [1,sqrt(3)]';

% Define AGA algorithm functions
unifun = @(x,f) deal(x,f); % Discard identical individuals (unimplemented)
fitfun = @(x) [norm(x-xa), norm(x-xb), norm(x-xc)]; % Fitness function (to be minimized)
mutfun = @(x,f) x + ranrange(-0.5,0.5,2); % Mutation: small random movement
repfun = @(x,y,fx,fy) (x+y)/2; % Reproduction: average of two individuals
ranfun = @() ranrange(0,10,2); % Random individual
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repeatability in the heuristic

% Execute Genetic Algorithm (GA)
[ bestInd, bestFit, nite, lastPop, lastFit, history ] = agap ( ...
    opts, np, goal, ng, N, npf, unifun, fitfun, ...
    mutfun, repfun, ranfun, prifun );

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('AGA \t\t%1.6f,%1.6f \t\t%1.6E,%1.6E,%1.6E\n',bestInd,bestFit);

