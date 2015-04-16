%% Example APS
% Find minima of a function with Particle Swarm (PS) algorithm
%
% Programmers:   Manel Soria         (UPC/ETSEIAT)
%                David de la Torre   (UPC/ETSEIAT)
%                Arnau Miro          (UPC/ETSEIAT)
% Date:          16/04/2015
% Revision:      1

%% APS

% Our test is a R^2->R function based on Rastrigin function.
% It is challenging because it has infinite local extrema, located at
% integer numbers (ie, 8,-9)
% The global minimum is at (1,1), and its value is 0
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define PS function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 0; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})

% Define PS parameters
nitemax = 50; % Maximum number of iterations
np = 50; % Population size
c1 = 10; % Local learning factor
c2 = 20; % Global learning factor
P1 = 1E-2; % Maximum step size (velocity) allowed in one iteration
goal = 1E-5; % Target fitness value

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define PS functions
fitfun = @(x) ras(x(1),x(2)); % Fitness function - TO BE MINIMIZED
ranfun = @() ranrange(-5,5,2); % Random individual
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repeatability in the PS

% Generate initial population
for i=1:np
    pop{i} = ranfun(); %#ok
    v{i} = 0.1 * ranfun(); %#ok
end;

% Execute Particle Swarm
[ bestIndAPS, bestFitAPS, nite, lastPopAPS, lastFitAPS, history ] = ...
    aps ( opts, pop, v, c1, c2, P1, nitemax, goal, ... 
    fitfun, prifun );

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1E-8,'Display','none');
[bestIndFMS,bestFitFMS] = fminsearch(fitfun,bestIndAPS,options);

% Display results of aps and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('APS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndAPS,bestFitAPS);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndFMS,bestFitFMS);

%% Fitness plot

% Get fitness history
if opts.nhist>1 && iscell(history) % Full history; get fitness values
    fithist = zeros(length(history),1);
    for i=1:length(history)
        fithist(i) = history{i,2}(1);
    end;
else fithist = history; % Simple history
end;

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

end;

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
    for iter=1:length(history)

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

        end;

        % Do events
        axis([-5 5 -5 5]);
        drawnow;
        
        % Wait
        pause(0.1);

        % Delete individuals
        if iter~=length(history) % Keep last frame
            for i=1:np, delete(ph{i}); end;
        end;

    end;

end;

