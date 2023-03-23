%% Example AGA
% Find minima of a function with Genetic Algorithm (GA)
%
% Programmers:   Manel Soria         (UPC/ETSEIAT)
%                David de la Torre   (UPC/ETSEIAT)
%                Arnau Miro          (UPC/ETSEIAT)
% Date:          23/11/2016
% Revision:      3
close all; clc;

%% AGA

% Define heuristic function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 1; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})

% Define AGA algorithm parameters
goal = -1E99; % Target fitness value
ng = 10; % Number of generations
np = 5000; % Population size
N = [1,... % Number of elites
    floor(np*0.1),... % Number of mutants
    floor(np*0.05),...% Number of newcomers
    floor(np*0.2)]; % Number of parents

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a and b
perfun = @(x) 2 * (x(1) + x(2)); % Perimeter

xa = [0,0]';
xb = [2,0]';
xc = [1,sqrt(3)]';
% xc = [1,0]';
% xc = [0,1]';
% xd = [2,1]';

xa = [0.5,0]';
xb = [1.5,0]';
xc = [1,0.5*sqrt(3)]';

% Define AGA algorithm functions
unifun = @(x,f) deal(x,f); % Discard identical individuals (unimplemented)
fitfun = @(x) [norm(x-xa), norm(x-xb), norm(x-xc)]; % Fitness function (to be minimized)
mutfun = @(x,f) x + ranrange(-0.1,0.1,2); % Mutation: small random movement
repfun = @(x,y,fx,fy) (x+y)/2; % Reproduction: average of two individuals
ranfun = @() ranrange(0,2,2); % Random individual
prifun = @(x) fprintf('%f %f ',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repeatability in the heuristic

% Execute Genetic Algorithm (GA)
[ bestInd, bestFit, nite, lastPop, lastFit, history ] = agap ( ...
    opts, np, goal, ng, N, unifun, fitfun, ...
    mutfun, repfun, ranfun, prifun );

% Display results of aga and fminsearch algorithms
fprintf('\nAlgorithm \tBest individual (x,y) \tValue\n');
fprintf('AGA \t\t%1.6f,%1.6f \t\t%1.6E\n',bestInd,bestFit);

%% Fitness plot

% Get fitness history
if opts.nhist>1 && iscell(history) % Full history; get fitness values
    fithist = zeros(length(history),3);
    for i=1:length(history)
        fithist(i,:) = history{i,2}(1,:);
    end
else % Simple history
    fithist = history;
end

% Plot fitness history
if ~isempty(fithist)

    % Create figure
    fh1 = figure('Position',[400,200,900,600]);

    % Beautify plot
    grid minor;
    title('Genetic Algorithm optimization | Rastrigin function');
    xlabel('Generation [#]');
    
    % Plot history
    yyaxis left; 
    ylabel('Best fitness function value #1 [log]');
    semilogy(fithist(:,1),'o-');

    % Plot history
    yyaxis right; 
    ylabel('Best fitness function value #2 [log]');
    semilogy(fithist(:,2),'x-');

end

%% Generations plot

% Only show generations when outputting full history
if opts.nhist>1 && iscell(history)

    % Create figure
    fh2 = figure('Position',[400,200,900,600]);
    view(0,90);
    hold on; 
    grid on;
    box on;
    axis equal;
    
    % Limits
    xlim([0,2]);
    ylim([0,2]);
    
    % Labels
    xlabel('Rectangle base [length]');
    ylabel('Rectangle height [length]');

    % Plot targets
    plot(xa(1),xa(2),'ro');
    plot(xb(1),xb(2),'ro');
    plot(xc(1),xc(2),'ro');
%     plot(xd(1),xd(2),'ro');
    
    % Colormap
    cm = jet(50);
    colormap(cm);
    hc = colorbar;
    title(hc, '# pareto front');

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
            if i<=ne, marker = 'v'; % Elites (rv)
            elseif i<=ne+nm, marker = 'o'; % Mutants (mo)
            elseif i<=ne+nm+nd, marker = 'x'; % Descendants (bx)
            else, marker = 's'; % Newcomers (ks)
            end
            
            % Select color
            idx_front = history{g,3}(i);
            if isnan(idx_front), c = 'k';
            else, c = cm(idx_front,:); end

            % Plot individual
            x = history{g,1}{i}(1);
            y = history{g,1}{i}(2);
            z = idx_front;
            ph{i} = plot3(x,y,z,marker,'MarkerSize',4,'color',c);

            % Save legend ticks
            if i==ne, lh(1) = ph{i}; % Elite
            elseif i==ne+nm, lh(2) = ph{i}; % Mutant
            elseif i==ne+nm+nd, lh(3) = ph{i}; % Descendant
            elseif i==ne+nm+nd+1, lh(4) = ph{i}; % Newcomer
            end

        end

        % Legend
        legend(lh(1:4),'Elites','Mutants','Descendants','Newcomers',...
            'Location','NorthEastOutside');

        % Do events
        drawnow;
        
        % Wait
        pause(1);

        % Delete individuals
        if g~=length(history) % Keep last frame
            for i=1:np, delete(ph{i}); end
        end

    end

end

