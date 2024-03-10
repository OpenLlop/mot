%% Example ADE
% Find minima of a function with Differential Evolution (DE)
%
% Programmers:   David de la Torre   (UPC/ETSEIAT)
%                Manel Soria         (UPC/ETSEIAT)
%                Arnau Miro          (UPC/ETSEIAT)
% Date:          23/11/2016
% Revision:      3

%% ADE

% Our test is a R^2->R function based on Rastrigin function.
% It is challenging because it has infinite local extrema, located at
% integer numbers (ie, 8,-9)
% The global minimum is at (1,1), and its value is 0
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define heuristic function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 1; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})

% Define ADE parameters
ng = 50; % Number of generations
np = 200; % Population size
N = [3,... % Number of elites
    floor(np*0.7)]; % Number of mutants
F = 0.1; % Mutation scaling factor
ms = 1; % Mutation strategy (see ade.m)
goal = 1E-5; % Target fitness value

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define ADE functions
unifun = @(x,f) deal(x,f); % Discard identical individuals (unimplemented)
fitfun = @(x) ras(x(1),x(2)); % Fitness function - TO BE MINIMIZED
mutfun = @(F,a,b,c) a + F * rand() * (b - c); % Mutation: random vector movement
ranfun = @() ranrange(-5,5,2); % Random individual
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repeatability in the heuristic algorithm

% Execute Differential Evolution (DE)
[ bestInd, bestFit, nite, lastPop, lastFit, history ] = ade ( ...
    opts, np, goal, ng, N, F, ms, unifun, fitfun, mutfun, ranfun, prifun );

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1E-8,'Display','none');
[bestIndFMS,bestFitFMS] = fminsearch(fitfun,bestInd,options);

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('ADE \t\t%1.6f,%1.6f \t\t%1.6E\n',bestInd,bestFit);
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

% Plot fitness histroy
if ~isempty(fithist)

    % Create figure
    fh1 = figure('Position',[400,200,900,600]);

    % Plot history
    semilogy(fithist,'o-');

    % Beautify plot
    grid minor;
    title('Differential Evolution optimization | Rastrigin function');
    xlabel('Generation [#]');
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

    % Population size
    ne = N(1); % Number of elites
    nm = N(2); % Number of mutants

    % Plot generations
    ph = cell(np,1); % Handles
    for g=1:ngens

        % Title
        title({['Differential Evolution optimization',...
            ' | Rastrigin function'];sprintf('Generation %03.0f',g)});
        
        % Plot individuals
        for i=1:np

            % Select plotting marker
            if i<=ne, marker = 'rv'; % Elites
            elseif i<=ne+nm, marker = 'mo'; % Mutants
            else, marker = 'ks'; % Newcomers
            end

            % Plot individual
            x = history{g,1}{i}(1);
            y = history{g,1}{i}(2);
            z = 100;
            ph{i} = plot3(x,y,z,marker,'MarkerSize',4);

            % Save legend ticks
            if i==ne, lh(1) = ph{i}; % Elite
            elseif i==ne+nm, lh(2) = ph{i}; % Mutant
            elseif i==ne+nm+1, lh(3) = ph{i}; % Newcomer
            end

        end

        % Legend
        legend(lh(1:3),'Elites','Mutants','Newcomers',...
            'Location','NorthEastOutside');

        % Do events
        drawnow;
        
        % Wait
        pause(1);

        % Delete individuals
        if g~=ngens % Keep last frame
            for i=1:np, delete(ph{i}); end
        end

    end

end

