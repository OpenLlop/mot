% Example AGA
% Find minima of a function with Genetic Algorithm
% Manel Soria, David de la Torre and Arnau Miro - ETSEIAT

% Clean-up
close all;
clear;

%% AGA

% Our test is a $R^2 \rightarrow R$ function based on Rastrigin function.
% It is challenging because it has infinite local extrema, located at
% integer numbers (ie, 8,-9)
% The global minimum is at (1,1), and its value is 0
ras = @(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Define GA function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 1; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})
opts.plotf = 1; % Plot fitness (0=none, 1=plot, 2=plot+save)
opts.plotp = 1; % Plot population (0=none, 1=plot, 2=plot+save)

% Define GA parameters
ng = 50; % Number of generations
np = 200; % Population size
N = [3,... % Number of elites
    floor(np*0.1),... % Number of mutants
    floor(np*0.05),...% Number of newcomers
    floor(np*0.2)]; % Number of parents
goal = 1E-5; % Target fitness value

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define GA functions
funique = @(x) x; % Discard identical individuals: currently not in use
fitfun = @(x) ras(x(1),x(2)); % Fitness function - TO BE MINIMIZED
mutfun = @(x,f) x + ranrange(-0.1,0.1,2); % Mutation: small random mov
repfun = @(x,y,fx,fy) (x+y)/2; % Reproduction: average
ranfun = @() ranrange(-5,5,2); % Random individual
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repeatability in the GA

% Execute Genetic Algorithm
[lastPopGA,bestFitGA,nite,history] = aga(opts,np,ng,N,goal,...
    funique,fitfun,mutfun,repfun,ranfun,prifun);

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1E-8,'Display','none');
[bestIndFMS,bestFitFMS] = fminsearch(fitfun,lastPopGA{1},options);

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('AGA \t\t%1.6f,%1.6f \t\t%1.6E\n',lastPopGA{1},bestFitGA);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6E\n',bestIndFMS,bestFitFMS);

%% Fitness plot

% Get fitness history
if opts.nhist>1 && iscell(history) % Full history; get fitness values
    histFitness = zeros(length(history),1);
    for i=1:length(history)
        histFitness(i) = history{i,2}(1);
    end;
else histFitness = history; % Simple history
end;

% Plot data
if ~isempty(histFitness)

    % Create figure
    fh1 = figure('Position',[400,200,900,600]);

    % Plot history
    semilogy(histFitness,'o-');

    % Beautify plot
    grid minor;
    title('Genetic Algorithm optimization | Rastrigin function');
    xlabel('Generation [#]');
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

    % Population size
    ne = N(1); % Number of elites
    nm = N(2); % Number of mutants
    nd = np - N(1) - N(2) - N(3); % Number of descendants

    % Plot generations
    ph = cell(np,1); % Handles
    for g=1:length(history)

        % Title
        title({'Genetic Algorithm optimization | Rastrigin function';...
            sprintf('Generation %03.0f',g)});
        
        % Plot individuals
        for i=1:np

            % Select plotting marker
            if i<=ne, marker = 'rv'; % Elites
            elseif i<=ne+nm, marker = 'mo'; % Mutants
            elseif i<=ne+nm+nd, marker = 'bx'; % Descendants
            else marker = 'ks'; % Newcomers
            end;

            % Plot individual
            x = history{g,1}{i}(1);
            y = history{g,1}{i}(2);
            z = 100;
            ph{i} = plot3(x,y,z,marker,'MarkerSize',4);

            % Save legend ticks
            if i==ne, lh(1) = ph{i}; % Elite
            elseif i==ne+nm, lh(2) = ph{i}; % Mutant
            elseif i==ne+nm+nd, lh(3) = ph{i}; % Descendant
            elseif i==ne+nm+nd+1, lh(4) = ph{i}; % Newcomer
            end;

        end;

        % Legend
        legend(lh,'Elites','Mutants','Descendants','Newcomers',...
            'Location','NorthEastOutside');

        % Do events
        drawnow;
        
        % Wait
        pause(1);

        % Delete individuals
        if g~=length(history) % Keep last frame
            for i=1:np, delete(ph{i}); end;
        end;

    end;

end;

