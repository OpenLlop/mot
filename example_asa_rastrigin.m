% Example ASA
% Find minima of a function with Simulated Annealing (SA) algorithm
% Manel Soria, David de la Torre and Arnau Miro - ETSEIAT

% Clean-up
close all;
clear;

%% ASA

% Our test is a R^2->R function based on Rastrigin function.
% It is challenging because it has infinite local extrema, located at
% integer numbers (ie, 8,-9)
% The global minimum is at (1,1), and its value is 0
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define SA function options
opts.ninfo = 10; % Verbosity level (print every # iterations)
opts.einfo = 0; % Print extended information
opts.label = 10; % Label (identification purposes)
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all data)

% Define SA parameters
nitemax = 50; % Maximum number of iterations
mu = 0.2; % Thermal transition probability parameter
goal = 1E-5; % Target fitness value

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define SA functions
fitfun = @(x) ras(x(1),x(2)); % Fitness function - TO BE MINIMIZED
mutfun = @(x,f) x + ranrange(-0.1,0.1,2); % Mutation: small random mov
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Initial guess
A0 = [2*rand(); 2*rand()];

% Execute Simulated Annealing
[ bestIndASA, bestFitASA, nite, lastPopASA, lastFitASA, history ] = ...
    asa ( opts, A0, nitemax, mu, goal, ...
    fitfun, mutfun, prifun );

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1e-8,'Display','none');
[bestIndFMS,bestFitFMS] = fminsearch(fitfun,bestIndASA,options);

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('ASA \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndASA,bestFitASA);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndFMS,bestFitFMS);

%% Plot fitness

% Get fitness history
if opts.nhist>1 % Full history; get fitness values
    fithist = zeros(length(history),1);
    for i=1:length(history)
        fithist(i) = history{i,6};
    end;
else fithist = history; % Simple history
end;

% Create figure
fh1 = figure('Position',[400,200,900,600]);
hold on;

% Plot history
semilogy(fithist,'o-');

% Beautify plot
grid minor;
title('Rastrigin function | Simulated Annealing optimization');
xlabel('Iteration [#]');
ylabel('Best fitness function value [log]');
hold off;


%% Plot iterations

% Only show generations when outputting full history
if opts.nhist>1 && iscell(history)

    % Create figure
    fh2 = figure('Position',[400,200,900,600]);

    % Plot rastrigin function
    [x,y] = meshgrid(-3:0.05:3,-3:0.05:3); z = ras(x,y);
    bh = surf(x,y,z,'LineStyle','none');
    colorbar('Location','EastOutside');
    view(0,90); hold on;
    
    % Virtual position for population
    z = 100;
    
    % Legend
    lh = plot3(0,0,-z,'rx',0,0,-z,'mo',0,0,-z,'co');
    legend(lh,'Best','A','B','Location','NorthEastOutside');

    % Plot generations
    ph = cell(length(history),3); % Handles
    for iter=1:length(history)

        % Title
        title({'Genetic Algorithm optimization | Rastrigin function';...
            sprintf('Iteration %03.0f',iter)});
        
        % Plot best
        x = history{iter,5}(1); y = history{iter,5}(2);
        ph{iter,1} = plot3(x,y,z,'rx','MarkerSize',8);

        % Plot A
        x = history{iter,1}(1); y = history{iter,1}(2);
        ph{iter,2} = plot3(x,y,z,'mo','MarkerSize',8);

        % Plot B
        x = history{iter,2}(1); y = history{iter,2}(2);
        ph{iter,3} = plot3(x,y,z,'co','MarkerSize',8);

        % Wait
        pause(0.1);

        % Delete individuals
        if iter~=length(history) % Keep last frame
            for i=1:3, delete(ph{iter,i}); end;
        end;

    end;

end;

