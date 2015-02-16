function [ lastpop, lastfit, nite, history ] = aga ( ninfo, label, ...
    pop, ng, N, goal, ...
    funique, fitfun, mutfun, repfun, ranfun, prifun )  

% 1-parfor que sigui opcional
% 2-si alguna funcio es empty que no la cridi si no es
% imprescindible (ie, print, mutacio)
% 3-abans de calcular, mirar si ja hem calculat
%
%-fem servir una funcio isequivalent(a,b) de l'usuari
% si torna 1 -> no es recalcula
% si torna 0 -> si es recalcula
% -si isequivalent es empty, internament fa servir isequuak
% va construint una llista de individuos coneguts, fins arribar a NCACHE
% si NCACHE=0, no fa res d'aixo (en algun cas sera lo millor ja que el 
% cost de buscarlo es creixent amn NCACHE*NP
                       
% Iterates to find minimum of a function using Genetic Algorithm
% (c) 2013 - Manel Soria - ETSEIAT - v1.01
% (c) 2015 - Manel Soria, David de la Torre - ETSEIAT - v1.02
%
% ninfo:    iteration control; prints every ninfo iterations 
% label:    integer number that precedes the prints in case output is to be
%           filtered
% pop:      list with initial population elements
% ng:       number of generations
% N:        population control parameters
%   N(1)    ne: number of elite individuals that remain unchanged
%   N(2)    nm: number of mutants
%   N(3)    nn: number of newcomers
%               the rest are descendants
%   N(4)    na: number of parents. The descendants are choosen 
%                   among the na best individuals 
%           The rest of individuals (ie: nn=length(pop)-ne+nm+nd) are
%                   newcommers, randomly choosen
% goal:     If function value is below goal, iterations are stopped
% 
% If there are less than nm-1 non-identical indivials, population is 
% considered degenerate and iterations stop
%
% Call back functions to be provided by user:
% funique:  Deletes repeated individuals in a population
%           Receives a population and returns a population
%           (a population is a list of individuals)
% fitfun:   Fitness function, given one individual returns its fitness
%           (RECALL that in this GA algoritm fitness is MINIMIZED)
% mutfun:   Mutation funcion, given one individual and its fitness,
%           mutfun should return a mutant individual. Fitness is given
%           in case mutation intensity is to be decreased when close to
%           the goal
% repfun:   Given two individuals, returns a descendant 
% ranfun:   Returns a random individual
% prifun:   Prints individual
%
% aga returns:
% lastpop:  list with the latest population sorted by fitness
% lastfit:  minimum value of fitfun found (from latest population)
% nite:     number of iterations performed
% history:  vector with the best value found after each iteration
       
% Build population if required
if isnumeric(pop) 
    NI = pop;
    pop = cell(1,NI);
    for i=1:NI
        pop{i} = ranfun(); % Create random individual
    end;
end;

% Population size
ne = N(1); % Number of elites
nm = N(2); % Number of mutants
nn = N(3); % Number of newcomers
na = N(4); % Number of breeders (selected from the best individuals)
ps = length(pop); % Population size
nd = ps - N(1) - N(2) - N(3); % Number of descendants

% History
history = [];

 % Safety checks
if nn<0
    error('aga: nn must be positive');
end;

% Iterate through generations
for g=1:ng
    
    % Save current generation index
    nite = g;

    % Clean population: remove repeated individuals
    pop = funique(pop); % Return unique individuals
    cps = length(pop); % Length of clean population

    % Avoid population degeneration (i.e., poor genetic pool)
    if cps<na % Clean population size is less than breeders size
        if info>0
            fprintf('GA label=%d degenerate population\n',label);
        end;
        break; % Break execution
    end

    % Repopulation: fill clean population pool with new individuals
    for i=cps+1:ps % Fill up to initial population size
        pop{i} = ranfun(); % Create new random individual
    end

    % Evaluate fitness function
    parfor i=1:ps
        fi(i) = feval(fitfun,pop{i});
    end;

    % Sort population individuals by their fitness level
    [fi,i] = sort(fi); % Sort fitness by increasing value (lower is best)
    pop = pop(i); % Sort population by their fitness value

    % Save fitness history
    history(g) = fi(1); %#ok

    % Show info if required
    if ninfo>1
        fprintf('GA label=%d g=%3d ng=%d best=%e ',label,g,ng,fi(1));
        if ~isempty(prifun)
            prifun(pop{1}); % Print best individual
        end;
        fprintf('\n');
    end;

    % Simulation end: either reached target fitness or max generations 
    if fi(1)<=goal || g>=ng
        
        % Save last iteration data
        lastpop = pop; % Save population
        lastfit = fi(1); % Save best fitness level
        
        % Show info if required
        if ninfo>0
            fprintf('GA label=%d best=%e ',label,lastfit);
            if ~isempty(prifun)
                prifun(pop{1}); % Print best individual
            end;
            if lastfit<goal % Goal achieved
                fprintf('goal=%e achieved !!\n',goal);
            else % Maximum generations reached (goal not achieved)
                fprintf('max. iterations reached, leaving\n');
            end;
        end;
        return; % Return from function
    end

    % Next generation:
    % <<[elites, mutants, descendants, newcomers]<<
    pop_next = cell(1,ne+nm+nd+nn); % Temp population
    k = 1; % Iteration index

    for i=1:ne % Elites
        pop_next{k} = pop{k}; % Copy
        k=k+1;
    end

    for i=1:nm % Mutants
        pop_next{k} = mutfun(pop{k},fi(k)); % Mutate
        k=k+1;
    end

    for i=1:nd % Descendants
        parentA = randi([1,na]); % Parent is choosen among np best 
        parentB = randi([1,na]); % Parent is choosen among np best 
        pop_next{k} = repfun(pop{parentA},pop{parentB}); % Breed
        k=k+1;
    end       

    for i=1:nn % Newcommers
        pop_next{k} = ranfun(); % Random individual
        k=k+1;
    end
    
    % Update population 
    pop = pop_next; % Update population
    clear('pop_next'); % Clear temp variable to conserve memory
    
end;

end

