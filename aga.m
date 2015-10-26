function [ bestind, bestfit, nite, lastpop, lastfit, history ] = ...
    aga( opts, pop, ng, N, goal, ...
    unifun, fitfun, mutfun, repfun, ranfun, prifun )
%AGA finds minimum of a function using Genetic Algorithm (GA)
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre   (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          14/04/2015
%Revision:      2
%
%Usage:         [bestind, bestfit, nite, lastpop, lastfit, history] = ...
%                   AGA( opts, pop, ng, N, goal, ...
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
%   unifun:     Deletes individuals considered to be too similar or equal in a population
%               Receives a population, sorted by fitness, and returns a population
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

% Get options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else ninfo = 1; end;
if isfield(opts,'label'), label = opts.label; else label = 0; end;
if isfield(opts,'dopar'), dopar = opts.dopar; else dopar = 0; end;
if isfield(opts,'nhist'), nhist = opts.nhist; else nhist = 1; end;

% Create history array
history = [];

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
if ne<0, error('aga: ne (number of elites) must be positive or zero'); end; 
if nm<0, error('aga: nm (number of mutants) must be positive or zero'); end;
if nn<0, error('aga: nn (number of newcomers) must be positive or zero'); end;
if na<0, error('aga: na (number of breeders) must be positive or zero'); end;
if nd<0, error('aga: nd (number of descendants) must be positive or zero'); end; 

if ninfo>0
    fprintf('GA label=%d np=%d population size\n',label,np);
    fprintf('GA label=%d ne=%d elite\n',label,ne);
    fprintf('GA label=%d nm=%d mutants\n',label,nm);
    fprintf('GA label=%d nn=%d newcommers\n',label,nn);
    fprintf('GA label=%d nd=%d descendants, from na=%d ancestors \n',label,nd,na);
end

% Iterate through generations
for g=1:ng

    % Evaluate fitness function
    fi = zeros(np,1); % Preallocate/clear fitness array    
    if dopar % Parallel execution
        parfor i=1:np, fi(i) = feval(fitfun,pop{i}); end;
    else % Serial execution
        for i=1:np, fi(i) = feval(fitfun,pop{i}); end;
    end;

    % Sort population individuals by their fitness level
    [fi,i] = sort(fi); % Sort fitness by increasing value (lower is best)
    pop = pop(i); % Sort population by fitness

    % Clean population: call the function that removes individuals too
    % similar (if it has been implemented)
    if ~isempty(unifun)
        [pop,diversity] = unifun(pop,fi); % Delete unique individuals, get a measure of population diversity
    else
        diversity=[];
    end
    ncleanpop = length(pop); % Length of clean population

    nnu=np-length(pop); % number of non unique individuals 
    
    % Avoid population degeneration (i.e., poor genetic pool)
    if ncleanpop<na % Stop if unique size is less than breeders size        
        % Show info
        if ninfo>0
            fprintf('GA label=%d degenerate population nnu=%d\n',label,nnu);
        end;
        
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


    % Save history
    if nhist>1 % Save full history {population,fitness}
        history{g,1} = pop; %#ok
        history{g,2} = fi; %#ok
    elseif nhist>0 % Save best fitness only
        history(g) = fi(1); %#ok
    end;
    
    % Show extended info
    if ninfo>1
        fprintf('GA label=%d g=%3d ng=%d nnu=%d best=%e ',label,g,ng,nnu,fi(1) ); 
        if ~isempty(diversity)
            fprintf('diver=%e',diversity); 
        end
        if ~isempty(prifun)
            fprintf(':');
            prifun(pop{1}); % Print best individual
        end;
        fprintf('\n');
    end;

    % Check if reached target fitness or max generations 
    if fi(1)<=goal || g>=ng % Target achieved
        
        % Save last iteration data
        bestind = pop{1}; % Save best individual
        bestfit = fi(1); % Save fitness level of last best individual
        nite = g; % Save current generation index
        lastpop = pop; % Save last population
        lastfit = fi; % Save last fitness values
        
        % Show info
        if ninfo>0
            fprintf('GA label=%d best=%e diver=%f bst:',label,bestfit,diversity);
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
    nextpop = cell(1,ne+nm+nd+nn); % Preallocate/clear next population
    k = 1; % Start individual counter index

    for i=1:ne % Elites
        nextpop{k} = pop{k}; % Copy elite into next generation
        k = k + 1; % Next individual
    end;

    for i=1:nm % Mutants
        if isempty(mutfun)
            nextpop{k} = pop{k}; % Do not mutate
        else
            nextpop{k} = mutfun(pop{k},fi(k)); % Mutate
        end;
        k = k + 1; % Next individual
    end;

    for i=1:nd % Descendants
        parentA = randi([1,na]); % Parent A is choosen among the na best 
        parentB = randi([1,na]); % Parent B is choosen among the na best 
        nextpop{k} = repfun(pop{parentA}, pop{parentB}, ...
            fi(parentA), fi(parentB)); % Breed individuals A and B
        k = k + 1; % Next individual
    end;

    for i=1:nn % Newcommers
        nextpop{k} = ranfun(); % Random individual
        k = k + 1; % Next individual
    end;
    
    % Update population 
    pop = nextpop;
    
end;

end

