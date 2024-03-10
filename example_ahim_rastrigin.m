%% Example Hybrid Islands Model + AGA/ADE/APS
% Find minima of a function with Hybrid Islands Model
%
% Programmers:   Manel Soria         (UPC/ETSEIAT)
%                David de la Torre   (UPC/ETSEIAT)
%                Arnau Miro          (UPC/ETSEIAT)
% Date:          29/12/2016
% Revision:      2

%% AHIM

% Our test is a R^2->R function based on Rastrigin function.
% It is challenging because it has infinite local extrema, located at
% integer numbers (ie, 8,-9).
% The global minimum is at (1,1), and its value is 0.
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define heuristic function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 1; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})

% Define Hybrid Islands Model parameters
ni = 3; % Number of islands
ngg = 6; % Number of global iterations
pops = [ni,40,40,40]; % Islands, and Population size of each island
nemi = 5; % Number of emigrants (must be lower than population of each island)
goal = -Inf; % Target fitness value
heufun = {@aga,@ade,@aps}; % Heuristic functions to use with the Hybrid Islands Model

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a and b

% Define optimization functions
unifun = @(x,f) deal(x,f); % Discard identical individuals: currently not in use
fitfun = @(x) ras(x(1),x(2)); % Fitness function (generic) (to be minimized)
mutfun = @(x,f) x + ranrange(-0.1,0.1,2); % Mutation for GA: small random movement
mutfun1 = @(F,a,b,c) a + F * rand() * (b - c); % Mutation for DE: small random movement
repfun = @(x,y,fx,fy) (x+y)/2; % Reproduction for GA: average
ranfun = @() ranrange(-5,5,2); % Random individual (generic)
rvlfun = @(vfact) vfact * ranrange(-5,5,2); % Random velocity for a PS individual
posfun = @(x,v) x + v; % Update PS particle position
velfun = @(v,x,xb,ib,c1,c2) v ... % Current PS velocity
    + c1*rand()*(xb - x) ... % Local PS learning factor
    + c2*rand()*(ib - x); % Global PS learning factor
vscfun = @(v,vmax) v*(norm(v)<vmax) ... % PS Velocity scaling not required
    + v*(vmax/norm(v))*(norm(v)>=vmax); % PS Limit excesive velocity
prifun = @(x) fprintf('%f %f',x(1),x(2)); % Print an individual (generic)

% Create DATA cell array
DATA = cell(ni,1);

% Assemble AGA data structure
DATA{1}{1} = 6; % Number of local generations
DATA{1}{2} = [3,... % Number of elites
    floor(pops(2)*0.1),... % Number of mutants
    floor(pops(2)*0.05),...% Number of newcomers
    floor(pops(2)*0.2)]; % Number of parents
DATA{1}{3} = unifun;
DATA{1}{4} = fitfun;
DATA{1}{5} = mutfun;
DATA{1}{6} = repfun;
DATA{1}{7} = ranfun;
DATA{1}{8} = prifun;

% Assemble ADE data structure
DATA{2}{1} = 6;
DATA{2}{2} = [3,... % Number of elites
    floor(pops(3)*0.7)]; % Number of mutants
DATA{2}{3} = 0.1; % Mutation scaling factor
DATA{2}{4} = 1; % Mutation strategy (see ade.m)
DATA{2}{5} = unifun;
DATA{2}{6} = fitfun;
DATA{2}{7} = mutfun1;
DATA{2}{8} = ranfun;
DATA{2}{9} = prifun;

% Assemble APS data structure
DATA{3}{1} = 50; % Maximum number of iterations
DATA{3}{2} = 0.2; % Step size (velocity) allowed in one iteration
DATA{3}{3} = 10; % Local learning factor
DATA{3}{4} = 20; % Global learning factor
DATA{3}{5} = 0.2; % Maximum step size (velocity) allowed in one iteration
DATA{3}{6} = fitfun;
DATA{3}{7} = posfun;
DATA{3}{8} = velfun;
DATA{3}{9} = vscfun;
DATA{3}{10} = rvlfun;
DATA{3}{11} = ranfun;
DATA{3}{12} = prifun;

% Randomize random seed
rng('shuffle'); % We don't want repeatability in the heuristic

% Execute Hybrid Islands Model
[ bestind, bestfit, nite, lastpop, lastfit, history ] = ahim ( ...
    opts, pops, ngg, nemi, goal, heufun, DATA );

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1e-8,'Display','none');
[bestIndFMS,bestFitFMS] = fminsearch(fitfun,bestind,options);

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('AGAI \t\t%1.6f,%1.6f \t\t%1.6E\n',bestind,bestfit);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndFMS,bestFitFMS);

%% Fitness plot

% Number of generations in history array
ngens = size(history,1);

% Get fitness history
if opts.nhist>1 && iscell(history) % Full history; get fitness values
    fithist = zeros(ngens,1);
    for g=1:ngens
        fithist(g) = history{g,ni+2};
    end
else % Simple history
    fithist = min(history,[],2);
end

% Plot data
if ~isempty(fithist)

    % Create figure
    fh1 = figure('Position',[400,200,900,600]);

    % Plot history
    semilogy(fithist,'o-');

    % Beautify plot
    grid minor;
    title(['Hybrid Islands Model optimization',...
        ' | Rastrigin function']);
    xlabel('Generation [#]');
    ylabel('Best fitness function value [log]');

end

