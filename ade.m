function [ bestind, bestfit, nite, lastpop, lastfit, history ] = ade ( ...
    opts, pop, goal, ng, N, F, ms, unifun, fitfun, mutfun, ranfun, prifun )
%ADE finds minimum of a function using Differential Evolution (DE)
%
%Programmers:   David de la Torre   (UPC/ETSEIAT)
%               Manel Soria         (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          10/05/2018
%Revision:      3
%
%Usage:         [bestind, bestfit, nite, lastpop, lastfit, history] = ...
%                   ADE ( opts, pop, goal, ng, N, F, ms, ...
%                   unifun, fitfun, mutfun, ranfun, prifun )
%
%Inputs:
%   opts:       function control parameters [struct] (optional)
%       ninfo:  verbosity level (0=none, 1=minimal, 2=extended)
%       label:  integer number that precedes the prints in case output is
%               to be filtered
%       dopar:  parallel execution of fitness function [1,0]
%       nhist:  save history (0=none, 1=fitness, 2=all{pop,fit})
%           0:  history = []
%           1:  history(ng) = bestfit(i)
%           2:  history{ng,1:2} = {pop,fitness}
%   pop:        list with initial population elements
%   goal:       if function value is below goal, iterations are stopped
%   ng:         maximum number of generations allowed
%   N:          population control parameters
%       N(1)        ne: number of elite individuals that remain unchanged
%       N(2)        nm: number of mutants
%                   The rest of individuals (ie: nn=length(pop)-ne+nm) are
%                   newcomers, randomly chosen
%   F:          mutation scaling factor (usually between 0-1)
%   ms:         mutation strategy
%       1:      DE/rand/1/exp -> p = mutfun(F,a,b,c)
%       2:      DE/rand/2/exp -> p = mutfun(F,a,b,c,d,e)
%       3:      DE/best/1/exp -> p = mutfun(F,best,b,c)
%       4:      DE/best/2/exp -> p = mutfun(F,best,b,c,d,e)
%       5:      DE/rand-to-best/1/exp -> p = mutfun(F,best,a,b,c)
%       6:      DE/rand-to-best/2/exp -> p = mutfun(F,best,a,b,c,d,e)
% 
%   If there are less than nm-1 non-identical individuals, population
%   is considered degenerate and iterations stop
%
%   Call back functions to be provided by the user:
%   unifun:     deletes repeated individuals in a population
%               receives a population+fitness and returns a
%               population+fitness (a population is a list of individuals)
%   fitfun:     fitness function, given one individual returns its fitness
%               (RECALL that in this DE algorithm fitness is MINIMIZED)
%   mutfun:     given X individuals and their fitnesses, returns a
%               combination of the individuals (X will depend on the
%               mutation strategy)
%   ranfun:     returns a random individual
%   prifun:     prints best individual
%
%Outputs:
%   bestind:    best individual from the last generation
%   bestfit:    fitness value of best individual from the last population
%   nite:       number of iterations (generations) performed
%   lastpop:    population of last generation
%   lastfit:    fitness values of last population
%   history:    array with saved history array

% Get options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else, ninfo = 1; end
if isfield(opts,'label'), label = opts.label; else, label = 0; end
if isfield(opts,'dopar'), dopar = opts.dopar; else, dopar = 0; end
if isfield(opts,'nhist'), nhist = opts.nhist; else, nhist = 1; end

% Create history array
history = [];

% Build population if required
if isnumeric(pop) % Population size is given as input
    np = pop; % Population size
    pop = cell(1,np); % Preallocate population variable
    for i=1:np % Fill population
        pop{i} = ranfun(); % Generate random individual
    end
end

% Population size
ne = N(1); % Number of elites
nm = N(2); % Number of mutants
np = length(pop); % Population size
nn = np - (ne + nm); % Number of newcomers

% Safety checks
if nm<=0, nm=1; end % Ensure there is at least one mutant
if nn<0, error('ADE number of elites/mutants exceeds population\n'); end

% Iterate through generations
for g=1:ng
    
    % Preallocate variables
    fi = zeros(np,1); % Preallocate/clear fitness array

    % Clean population: remove repeated individuals
    [pop,fi] = unifun(pop,fi); % Return unique individuals
    pop = pop(~cellfun('isempty',pop)); % Remove empty individuals
    fi = fi(~cellfun('isempty',pop)); % Remove empty fitness
    ncleanpop = length(pop); % Length of clean population

    % Avoid population degeneration (i.e., poor genetic pool)
    if ncleanpop<nm % Clean population size is less than mutants size
        
        % Show info
        if ninfo>0
            fprintf('ADE label=%d degenerate population, leaving\n',label);
        end
        
        % Save last iteration data
        bestind = pop{1}; % Save best individual
        bestfit = fi(1); % Save fitness level of last best individual
        nite = g; % Save current generation index
        lastpop = pop; % Save last population
        lastfit = fi; % Save last fitness values
        
        % Stop iterating
        break;
        
    end

    % Repopulation: fill clean population pool with new individuals
    for i=ncleanpop+1:np % Fill up to initial population size
        pop{i} = ranfun(); % Create new random individual
    end

    % Evaluate fitness function
    if dopar % Parallel execution
        parfor i=1:np, fi(i) = feval(fitfun,pop{i}); end
    else % Serial execution
        for i=1:np, fi(i) = feval(fitfun,pop{i}); end
    end

    % Sort population individuals by their fitness level
    [fi,i] = sort(fi); % Sort fitness by increasing value (lower is best)
    pop = pop(i); % Sort population by fitness

    % Save history
    if nhist>1 % Save full history {population,fitness}
        history{g,1} = pop; %#ok
        history{g,2} = fi; %#ok
    elseif nhist>0 % Save best fitness only
        history(g) = fi(1); %#ok
    end
    
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
            fprintf('ADE label=%d nite=%2d fitbest=%f',label,nite,bestfit);
            if ~isempty(prifun), fprintf(' best= '); prifun(bestind); end
            if bestfit<goal % Goal achieved
                fprintf(' goal=%e achieved, leaving\n',goal);
            else % Maximum generations reached (goal not achieved)
                fprintf(' max. iterations reached, leaving\n');
            end
        end
        
        % Stop iterating
        break;
        
    end

    % Show extended info
    if ninfo>1
        fprintf('ADE label=%d g=%2d fitbest=%e',label,g,fi(1));
        if ~isempty(prifun), fprintf(' best='); prifun(pop{1}); end
        fprintf('\n');
    end
    
    % Compute population for next generation:
    % <<[elites, mutants, newcomers]<<
    nextpop = cell(1,ne+nm+nn); % Preallocate/clear next population
    k = 1; % Start individual counter index

    for i=1:ne % Elites
        nextpop{k} = pop{k}; % Copy elite into next generation
        k = k + 1; % Next individual
    end

    for i=1:nm % Mutants
        
        % Select mutation strategy
        switch ms
            
            case 1 % DE/rand/1/exp

                % Individuals are chosen among the nm best
                p = num2cell(randperm(nm,3));
                [a,b,c] = p{:};

                % Mutate individual
                nextpop{k} = mutfun(F,pop{a},pop{b},pop{c});
                
            case 2 % DE/rand/1/exp

                % Individuals are chosen among the nm best
                p = num2cell(randperm(nm,5));
                [a,b,c,d,e] = p{:};

                % Mutate individual
                nextpop{k} = mutfun(F,pop{a},...
                    pop{b},pop{c},pop{d},pop{e});
                
            case 3 % DE/best/1/exp

                % Individuals are chosen among the nm best
                p = num2cell(randperm(nm,2));
                [b,c] = p{:};

                % Mutate individual
                nextpop{k} = mutfun(F,pop{1},pop{b},pop{c});
                
            case 4 % DE/best/2/exp

                % Individuals are chosen among the nm best
                p = num2cell(randperm(nm,4));
                [b,c,d,e] = p{:};

                % Mutate individual
                nextpop{k} = mutfun(F,pop{1},pop{b},...
                    pop{c},pop{d},pop{e});
                
            case 5 % DE/rand-to-best/1/exp

                % Individuals are chosen among the nm best
                p = num2cell(randperm(nm,3));
                [a,b,c] = p{:};

                % Mutate individual
                nextpop{k} = mutfun(F,pop{1},...
                    pop{a},pop{b},pop{c});
                
            case 6 % DE/rand-to-best/2/exp

                % Individuals are chosen among the nm best
                p = num2cell(randperm(nm,5));
                [a,b,c,d,e] = p{:};

                % Mutate individual
                nextpop{k} = mutfun(F,pop{1},pop{a},...
                    pop{b},pop{c},pop{d},pop{e});
                
            otherwise % DE/rand/1/exp
                
                % Individuals are chosen among the nm best
                p = num2cell(randperm(nm,3));
                [a,b,c] = p{:};

                % Mutate individual
                nextpop{k} = mutfun(F,pop{a},pop{b},pop{c});
                
        end
        
        % Evaluate fitness of tentative new individual
        nextfit = feval(fitfun,nextpop{i});
        
        % Decide if new individual is accepted into new population
        if nextfit > fi(k) % New individual has worse fitness; rejected
            nextpop{k} = pop{k}; % Recover old individual
        end
        
        % Next individual
        k = k + 1;
        
    end

    for i=1:nn % Newcommers
        nextpop{k} = ranfun(); % Random individual
        k = k + 1; % Next individual
    end
    
    % Update population
    pop = nextpop;
    
end

end

