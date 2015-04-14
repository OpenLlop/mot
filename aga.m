function [ bestind, bestfit, nite, lastpop, lastfit, history ] = ...
    aga ( opts, pop, ng, N, goal, ...
    unifun, fitfun, mutfun, repfun, ranfun, prifun )
%AGA Finds minimum of a function using Genetic Algorithm (GA)
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre	(UPC/ETSEIAT)
%Date:          14/04/2015
%Revision:      2
%
%Usage:         [bestind, bestfit, nite, lastpop, lastfit, history] = ...
%                   aga ( opts, pop, ng, N, goal, ...
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
%   pop:        list with initial population elements
%   ng:         maximum number of generations allowed
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
%	Call back functions to be provided by the user:
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
%   bestind:    best individual from the last generation
%   bestfit:    fitness value of best individual from the last population
%   nite:       number of iterations (generations) performed
%   lastpop:    population of last generation
%   lastfit:    fitness values of last population
%   history:    array with saved history array

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

% Declare history array, if required
if nhist>0 || nargout>3, history = []; end;

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
        bestind = pop{1}; % Save best individual
        bestfit = fi(1); % Save fitness level of last best individual
        nite = g; % Save current generation index
        lastpop = pop; % Save last population
        lastfit = fi; % Save last fitness values
        
        % Stop iterating
        break;
        
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
        bestind = pop{1}; % Save best individual
        bestfit = fi(1); % Save fitness level of last best individual
        nite = g; % Save current generation index
        lastpop = pop; % Save last population
        lastfit = fi; % Save last fitness values
        
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
        
        % Stop iterating
        break;
        
    end;

    % Compute population for next generation:
    % <<[elites, mutants, descendants, newcomers]<<
    nextpop = cell(1,ne+nm+nd+nn); % Temp population
    k = 1; % Iteration index

    for i=1:ne % Elites
        nextpop{k} = pop{k}; % Copy
        k = k + 1;
    end;

    for i=1:nm % Mutants
        if isempty(mutfun), nextpop{k} = pop{k}; % Do not mutate
        else nextpop{k} = mutfun(pop{k},fi(k)); % Mutate
        end;
        k = k + 1;
    end;

    for i=1:nd % Descendants
        parentA = randi([1,na]); % Parent is choosen among np best 
        parentB = randi([1,na]); % Parent is choosen among np best 
        nextpop{k} = repfun(pop{parentA}, pop{parentB}, ...
            fi(parentA), fi(parentB)); % Breed
        k = k + 1;
    end;

    for i=1:nn % Newcommers
        nextpop{k} = ranfun(); % Random individual
        k = k + 1;
    end;
    
    % Update population 
    pop = nextpop; % Update population
    clear('nextpop'); % Clear temp variable to conserve memory
    
end;

end

