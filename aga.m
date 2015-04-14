function [ lastpop, bestfit, nite, history ] = aga ( opts, ...
    pop, ng, N, goal, ...
    unifun, fitfun, mutfun, repfun, ranfun, prifun )
%AGA Finds minimum of a function using Genetic Algorithm (GA)
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre	(UPC/ETSEIAT)
%Date:          14/04/2015
%Revision:      2
%
%Usage:         [lastpop, bestfit, nite, history] = aga ( opts, ...
%                   pop, ng, N, goal, ...
%                   unifun, fitfun, mutfun, repfun, ranfun, prifun )
%
%Inputs:
%   opts:       function control parameters [struct] (optional)
%       ninfo:  verbosity level (0=none, 1=minimal, 2=extended)
%       label:  integer number that precedes the prints in case output is
%               to be filtered
%       dopar:  parallel execution of fitness function [1,0]
%       nhist:  save history (0=none, 1=fitness, 2=all{pop,fit})
%                   0: history = []
%                   1: history(ng) = bestfit(i)
%                   2: history{ng,1:2} = {pop,fitness}
%       plotf:  plot fitness (0=none, 1=plot, 2=plot+save)
%                   0: no plot
%                   1: plot on each generation
%                   2: plot on each generation + save plot to file
%       plotp:  plot population (0=none, 1=plot, 2=plot+save)
%                   0: no plot
%                   1: plot on each generation
%                   2: plot on each generation + save plot to file
%   pop:        list with initial population elements
%   ng:         number of generations
%   N:          population control parameters
%       N(1)    ne: number of elite individuals that remain unchanged
%       N(2)    nm: number of mutants
%       N(3)    nn: number of newcomers
%               The rest are descendants
%       N(4)    na: number of parents. The descendants are choosen among
%               the na best individuals
%               The rest of individuals (ie: nn=length(pop)-ne+nm+nd) are
%               newcommers, randomly choosen
%	goal:       If function value is below goal, iterations are stopped
% 
%	If there are less than nm-1 non-identical indivials, population is 
%	considered degenerate and iterations stop
%
%	Call back functions to be provided by user:
%   unifun:     Deletes repeated individuals in a population
%               Receives a population and returns a population
%               (a population is a list of individuals)
%	fitfun:     Fitness function, given one individual returns its fitness
%               (RECALL that in this GA algoritm fitness is MINIMIZED)
%	mutfun:     Mutation funcion, given one individual and its fitness,
%               mutfun should return a mutant individual. Fitness is given
%               in case mutation intensity is to be decreased when close
%               to the goal
%	repfun:     Given two individuals and their fitnesses, returns a
%               descendant
%   ranfun:     Returns a random individual
%   prifun:     Prints individual
%
%Outputs:
%   lastpop:    population from the last generation
%   bestfit:    fitness value of best individual from the last population
%   nite:       number of iterations (generations) performed
%   history:    array with the best fitness value found after each
%               iteration

% TODO:
% 2-si alguna funcio es empty que no la cridi si no es
% imprescindible (ie, print, mutacio)
% 3-abans de calcular, mirar si ja hem calculat (cache)
%
%-fem servir una funcio isequivalent(a,b) de l'usuari
% si torna 1 -> no es recalcula
% si torna 0 -> si es recalcula
% -si isequivalent es empty, internament fa servir isequuak
% va construint una llista de individuos coneguts, fins arribar a NCACHE
% si NCACHE=0, no fa res d'aixo (en algun cas sera lo millor ja que el 
% cost de buscarlo es creixent amn NCACHE*NP

% Set options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else ninfo = 1; end;
if isfield(opts,'label'), label = opts.label; else label = 0; end;
if isfield(opts,'dopar'), dopar = opts.dopar; else dopar = 0; end;
if isfield(opts,'nhist'), nhist = opts.nhist; else nhist = 1; end;
if isfield(opts,'plotf'), plotf = opts.plotf; else plotf = 0; end;
if isfield(opts,'plotp'), plotp = opts.plotp; else plotp = 0; end;

% Declare history variable, if required
if nhist>0, history = []; end;

% Build population if required
if isnumeric(pop) % Population size is given as input
    np = pop; % Population size
    pop = cell(1,np); % Preallocate population variable
    for i=1:np % Fill population
        pop{i} = ranfun(); % Generate random individual
    end;
end;

% Population size
ne = N(1); % Number of elites
nm = N(2); % Number of mutants
nn = N(3); % Number of newcomers
na = N(4); % Number of breeders (selected from the best individuals)
np = length(pop); % Population size
nd = np - N(1) - N(2) - N(3); % Number of descendants

% Safety checks
if nn<0, error('aga: nn (number of newcomers) must be positive'); end;

% Iterate through generations
for g=1:ng
    
    % Save current generation index
    nite = g;
    
    % Preallocate vars
    fi = zeros(np,1); % Preallocate var

    % Clean population: remove repeated individuals
    pop = unifun(pop); % Return unique individuals
    ncleanpop = length(pop); % Length of clean population

    % Avoid population degeneration (i.e., poor genetic pool)
    if ncleanpop<na % Clean population size is less than breeders size
        
        % Info
        if ninfo>0
            fprintf('GA label=%d degenerate population\n',label);
        end;
        
        % Save last iteration data
        lastpop = pop; % Save last population
        bestfit = fi(1); % Save fitness level of last best individual
        
        % Return from function
        return;
        
    end;

    % Repopulation: fill clean population pool with new individuals
    for i=ncleanpop+1:np % Fill up to initial population size
        pop{i} = ranfun(); % Create new random individual
    end;

    % Evaluate fitness function
    if dopar % Parallel execution
        parfor i=1:np, fi(i) = feval(fitfun,pop{i}); end;
    else % Serial execution
        for i=1:np, fi(i) = feval(fitfun,pop{i}); end;
    end;

    % Sort population individuals by their fitness level
    [fi,i] = sort(fi); % Sort fitness by increasing value (lower is best)
    pop = pop(i); % Sort population by their fitness value

    % Save history
    if nhist>1 % Save full history {population,fitness}
        history{g,1} = pop; %#ok
        history{g,2} = fi; %#ok
    elseif nhist>0 % Save best fitness only
        history(g) = fi(1); %#ok
    end;
    
    % Plot fitness if required
    if plotf>0
       if ~exist('fhfit','var'), fhfit = initFitFigure(); end;
       plotFitHistory(fhfit); % Plot fitness history
    end;

    % Plot generation if required
    if plotp>0 && length(pop{1})>1
       if ~exist('fhpop','var'), [fhpop, ph] = initGenFigure(); end;
       plotPopCurrent(fhpop); % Plot current generation
    end;

    % Show info if required
    if ninfo>1
        fprintf('GA label=%d g=%3d ng=%d best=%e ',label,g,ng,fi(1));
        if ~isempty(prifun)
            prifun(pop{1}); % Print best individual
        end;
        fprintf('\n');
    end;

    % Check if reached target fitness or max generations 
    if fi(1)<=goal || g>=ng % Target achieved; end simulation
        
        % Save last iteration data
        lastpop = pop; % Save last population
        bestfit = fi(1); % Save fitness level of last best individual
        
        % Show info if required
        if ninfo>0
            fprintf('GA label=%d best=%e ',label,bestfit);
            if ~isempty(prifun)
                prifun(pop{1}); % Print best individual
            end;
            if bestfit<goal % Goal achieved
                fprintf('goal=%e achieved !!\n',goal);
            else % Maximum generations reached (goal not achieved)
                fprintf('max. iterations reached, leaving\n');
            end;
        end;
        
        % Return from function
        return;
        
    end;

    % Compute population for next generation:
    % <<[elites, mutants, descendants, newcomers]<<
    nextpop = cell(1,ne+nm+nd+nn); % Temp population
    k = 1; % Iteration index

    for i=1:ne % Elites
        nextpop{k} = pop{k}; % Copy
        k=k+1;
    end;

    for i=1:nm % Mutants
        if isempty(mutfun), nextpop{k} = pop{k}; % Do not mutate
        else nextpop{k} = mutfun(pop{k},fi(k)); % Mutate
        end;
        k=k+1;
    end;

    for i=1:nd % Descendants
        parentA = randi([1,na]); % Parent is choosen among np best 
        parentB = randi([1,na]); % Parent is choosen among np best 
        nextpop{k} = repfun(pop{parentA}, pop{parentB}, ...
            fi(parentA), fi(parentB)); % Breed
        k=k+1;
    end;

    for i=1:nn % Newcommers
        nextpop{k} = ranfun(); % Random individual
        k=k+1;
    end;
    
    % Update population 
    pop = nextpop; % Update population
    clear('nextpop'); % Clear temp variable to conserve memory
    
end;


%% Auxiliary functions

    % Create fitness plot figure
    function [fh] = initFitFigure()

        % Figure Sizing
        SS = get(0,'ScreenSize'); % Get User's Screen Size
        figW = 800; % Figure Width
        figH = 600; % Figure Height
        figW0 = SS(3)/2 - figW/2; % Initial x-Coordinate
        figH0 = SS(4)/2 - figH/2; % Initial y-Coordinate

        % Create figure
        fh = figure('Position',[figW0,figH0,figW,figH],...
            'PaperSizeMode','auto');
        
        % Plot settings
        hold on; % Hold figure
        box on; grid minor;
        title('Genetic Algorithm optimization');
        xlabel('Generation [#]');
        ylabel('Best fitness function value');
        
    end

    % Plot fitness history
    function plotFitHistory(fh)
        
        % Set current figure
        figure(fh);
        
        % Plot history
        plot(g,fi(1),'b*-');
        
        % Do events
        drawnow;
        
        % Save plot to file
        if plotf>1
            if ~isdir('output'), mkdir('output'); end; % Create dir
            print(fh,'-dpng','-r300',fullfile('output',...
                sprintf('AGA_FIT_%d',label))); % Print plot
        end;

    end

    % Create population plot figure
    function [fh, ph] = initGenFigure()

        % Figure Sizing
        SS = get(0,'ScreenSize'); % Get User's Screen Size
        figW = 800; % Figure Width
        figH = 600; % Figure Height
        figW0 = SS(3)/2 - figW/2; % Initial x-Coordinate
        figH0 = SS(4)/2 - figH/2; % Initial y-Coordinate

        % Create figure
        fh = figure('Position',[figW0,figH0,figW,figH],...
            'PaperSizeMode','auto');
        
        % Plot settings
        hold on; % Hold figure
        view(0,90); box on; grid minor;

        % Create plot handles
        ph = cell(np,1);
        
    end

    % Plot current generation
    function plotPopCurrent(fh)

        % Set current figure
        figure(fh);
        
        % Title
        title({'Genetic Algorithm optimization';...
            sprintf('Generation %03d',g)});

        % Delete old individuals
        if g~=1, for ii=1:np, delete(ph{ii}); end; end;

        % Plot new individuals
        for ii=1:np

            % Select plotting marker (elite, mutant, normal, newcomer)
            if ii<=ne, marker = 'rv'; % Elites range
            elseif ii<=ne+nm, marker = 'mo'; % Mutants range
            elseif ii<=ne+nm+nd, marker = 'bx'; % Descendants range
            else marker = 'ks'; % Newcomers range
            end;

            % Plot individual
            x = pop{ii}(1); % X
            y = pop{ii}(2); % Y
            z = fi(ii); % Fitness
            ph{ii} = plot3(x,y,z,marker,'MarkerSize',4); % Plot 3D
            
            % Save legend ticks
            if ii==ne, lh(1) = ph{ii}; % Elite
            elseif ii==ne+nm, lh(2) = ph{ii}; % Mutant
            elseif ii==ne+nm+nd, lh(3) = ph{ii}; % Descendant
            elseif ii==ne+nm+nd+1, lh(4) = ph{ii}; % Newcomer
            end;

        end;
        
        % Legend
        legend(lh,'Elites','Mutants','Descendants','Newcomers',...
            'Location','NorthEastOutside');

        % Do events
        drawnow;
        
        % Save plot to file
        if plotp>1
            if ~isdir('output'), mkdir('output'); end; % Create dir
            print(fh,'-dpng','-r300',fullfile('output',...
                sprintf('AGA_POP_%d_%03d',label,g))); % Print plot
        end;

    end

end

