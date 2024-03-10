%% Example APS
% Find minima of a function with Particle Swarm (PS) algorithm
%
% Programmers:   Manel Soria         (UPC/ETSEIAT)
%                David de la Torre   (UPC/ETSEIAT)
%                Arnau Miro          (UPC/ETSEIAT)
% Date:          29/12/2016
% Revision:      4

%% APS

% Our test is a R^2->R function based on Rastrigin function.
% It is challenging because it has infinite local extrema, located at
% integer numbers (ie, 8,-9)
% The global minimum is at (1,1), and its value is 0
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define heuristic function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 0; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})

% Define APS parameters
nitemax = 50; % Maximum number of iterations
np = 50; % Population size
c1 = 10; % Local learning factor
c2 = 20; % Global learning factor
vmax = 0.2; % Maximum step size (velocity) allowed in one iteration
goal = 1E-5; % Target fitness value

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define APS functions
fitfun = @(x) ras(x(1),x(2)); % Fitness function (to be minimized)
ranfun = @() ranrange(-5,5,2); % Random individual
rvlfun = @() 0.1 * ranrange(-5,5,2); % Random individual
posfun = @(x,v) x + v; % Update particle position
velfun = @(v,x,xb,ib,c1,c2) v ... % Current velocity
    + c1*rand()*(xb - x) ... % Local learning factor
    + c2*rand()*(ib - x); % Global learning factor
vscfun = @(v,vmax) v*(norm(v)<vmax) ... % Velocity scaling not required
    + v*(vmax/norm(v))*(norm(v)>=vmax); % Limit excesive velocity
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repeatability in the heuristic

% Generate initial population
pop = cell(np,1);
v = cell(np,1);
for i=1:np
    pop{i} = ranfun();
    v{i} = rvlfun();
end

% Execute Particle Swarm
[ bestInd, bestFit, nite, lastPop, lastFit, history ] = aps ( ...
    opts, pop, goal, nitemax, v, c1, c2, vmax, fitfun, posfun, ...
    velfun, vscfun, rvlfun, ranfun, prifun );

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1E-8,'Display','none');
[bestIndFMS,bestFitFMS] = fminsearch(fitfun,bestInd,options);

% Display results of aps and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('APS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestInd,bestFit);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndFMS,bestFitFMS);

%% Fitness plot

% Number of generations in history array
ngens = size(history,1);

% Get fitness history
if opts.nhist>1 && iscell(history) % Full history; get fitness values
    fithist = zeros(ngens,1);
    for i=1:ngens
        fithist(i) = history{i,2}(1);
    end
else % Simple history
    fithist = history;
end

% Plot data
if ~isempty(fithist)

    % Create figure
    fh1 = figure('Position',[400,200,900,600]);

    % Plot history
    semilogy(fithist,'o-');

    % Beautify plot
    grid minor;
    title('Particle Swarm optimization | Rastrigin function');
    xlabel('Iteration [#]');
    ylabel('Best fitness function value [log]');

end

%% Generations plot

% Only show generations when outputting full history
if opts.nhist>1 && iscell(history)

    % Create figure
    fh2 = figure('Position',[400,200,900,600]);

    % Plot rastrigin function
    [x,y] = meshgrid(-5:0.05:5,-5:0.05:5); z = ras(x,y);
    bh = surf(x,y,z,'LineStyle','none');
    colorbar('Location','EastOutside');
    view(0,90); hold on;

    % Plot generations
    ph = cell(np,1); % Handles
    for iter=1:ngens

        % Title
        title({'Particle Swarm optimization | Rastrigin function';...
            sprintf('Iteration %03.0f',iter)});
        
        % Plot individuals
        for i=1:np

            % Plot individual
            x = history{iter,1}{i}(1);
            y = history{iter,1}{i}(2);
            z = 100;
            ph{i} = plot3(x,y,z,'ro','MarkerSize',6);

        end

        % Do events
        axis([-5 5 -5 5]);
        drawnow;
        
        % Wait
        pause(0.1);

        % Delete individuals
        if iter~=ngens % Keep last frame
            for i=1:np, delete(ph{i}); end
        end

    end

end

