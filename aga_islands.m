function [ bestind, bestfit, nite, lastpop, lastfit, history ] = ...
    aga_islands( opts, pops, ngg, nemi, ng, N, goal, ...
    unifun, fitfun, mutfun, repfun, ranfun, prifun )
%AGA_ISLANDS finds minimum of a function using Genetic Algorithm (GA) with
%Islands
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre   (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          14/04/2015
%Revision:      2
%
%Usage:         [ bestind, bestfit, nite, lastpop, lastfit, history ] = ...
%                   AGA_ISLANDS( opts, pops, ngg, nemi, ng, N, goal, ...
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
%                   1: history(ngg,ni) = bestfit(i)
%                   2: history{ngg,ni} = [hist,nite] (refer to aga)
%   pop:        list with initial population elements
%   ngg:        number of global generations
%   nemi:       number of emigrations between islands per global generation
%   ng:         number of generations for each island per global generation
%   N:          population control parameters
%       N(1)    ne: number of elite individuals that remain unchanged
%       N(2)    nm: number of mutants
%       N(3)    nn: number of newcomers
%               The rest are descendants
%       N(4)    na: number of parents. The descendants are choosen among 
%               the na best individuals
%               The rest of individuals (ie: nn=length(pop)-ne+nm+nd) are
%               newcommers, randomly choosen
%   goal:       If function value is below goal, iterations are stopped
% 
%   If there are less than nm-1 non-identical indivials, population is 
%   considered degenerate and iterations stop
%
%   Call back functions to be provided by the user:
%   unifun:     Deletes repeated individuals in a population
%               Receives a population and returns a population
%               (a population is a list of individuals)
%   fitfun:     Fitness function, given one individual returns its fitness
%               (RECALL that in this GA algoritm fitness is MINIMIZED)
%   mutfun:     Mutation funcion, given one individual and its fitness,
%               mutfun should return a mutant individual. Fitness is given
%               in case mutation intensity is to be decreased when close
%               to the goal
%   repfun:     Given two individuals and their fitnesses, returns a
%               descendant
%   ranfun:     Returns a random individual
%   prifun:     Prints individual
%
%Outputs:
%   bestind:    best individual (among all the islands)
%   bestfit:    fitness value of best individual
%   nite:       number of global iterations (generations) performed
%   lastpop:    list with last populations of each island
%   lastfit:    best fitness values of last population of each island
%   history:    array with saved global history array

% Get options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else ninfo = 0; end;
if isfield(opts,'label'), label = opts.label; else label = 0; end;
if isfield(opts,'dopar'), dopar = opts.dopar; else dopar = 0; end;
if isfield(opts,'nhist'), nhist = opts.nhist; else nhist = 0; end;

% Create history array
history = [];

% Build population
if (~iscell(pops)) % Only population size is given                                   
    ni = pops(1); % Number of islands
    np = pops(2); % Population of an island
    pops = cell(1,ni); % Preallocate var
    for island=1:ni % Fill islands with population
        for i=1:np % Fill population with individuals
            pops{island}{i} = ranfun(); % Create random individual
        end;
    end;
else % Initial population is given. Check structural consistency
    ss = size(pops); % Size (m,n) of input population
    if (ss(1)~=1) % Error in pops shape
        error('Global population shape must be: pops{1,ni}');
    end;
    ni = ss(2); % Number of islands
    for i=1:ni
        ss = size(pops{i}); % Size 
        if (ss(1)~=1), error('Population shape must be: pop{1,ni}'); end;
        if i==1, np = ss(2); % Population length
        else if np~=ss(2), error('Population length mismatch'); end;
        end;
    end;
end;

% Show info
if ninfo>0
    fprintf('aga_islands begin ni=%d np=%d ngg=%d ng=%d\n',ni,np,ngg,ng);
end;

% Iterate through generations
for gg=1:ngg;
    
    % Bestfit of each island
    bestfits = zeros(ni,1);
    
    % Evolve each island separately
    for island=1:ni
        
        % AGA island options
        iopts.ninfo = ninfo-1;
        iopts.label = label + gg*1000 + island;
        iopts.dopar = dopar;
        iopts.nhist = nhist;
        
        % Execute AGA (Genetic Algorithm)
        [ibestind, ibestfit, initer, ilastpop, ~, hist] = aga(iopts, ...
            pops{island}, ng, N, goal,...
            unifun, fitfun, mutfun, repfun, ranfun, prifun);
        
        % Save values at the end of local island iteration
        pops{island} = ilastpop; % Save last population
        bestfits(island) = ibestfit; % Save best fitness value

        % Save history
        if nhist>1 % Save full history for each island
            history{gg,island} = {hist,initer}; %#ok
        elseif nhist>0 % Save best fitness only
            history(gg,island) = ibestfit; %#ok
        end;

        % Show extended info
        if ninfo>1
            fprintf(['aga_islands label= %d end ite global= %d ',...
                'island= %d fbest= %e'],label,gg,island,bestfits(island));
            if ~isempty(prifun) % Print best individual
                fprintf(' best: '); prifun(ibestind);
            end;
            fprintf('\n');
        end;
        
    end;
    
    % Find best individual of all islands and island where it lives
    [fbest, ibest] = min(bestfits);
    
    % Show info
    if ninfo>0
        fprintf('aga_islands label= %d gg= %d fbest= %8.3e island= %d',...
            label,gg,fbest,ibest);
        if ~isempty(prifun)
            fprintf(' best: '); prifun(pops{ibest}{1});
        end;
        fprintf('\n');
    end;

    % Save history
    if nhist>1 % Save full history for each island
        history{gg,ni+1} = pops{ibest}{1}; %#ok % Best individual
        history{gg,ni+2} = fbest; %#ok % Fitness of best individual
    end;

    % Check if reached target fitness or max generations 
    if fbest<=goal || gg>=ngg % Target achieved; end simulation
        
        % Show info if required
        if ninfo>1
            fprintf('aga_islands stoping because ');
            if gg==ngg
                fprintf(['maximum number of global iterations %d ',...
                    'has been reached \n'], ng);
            else
                fprintf('goal %e has been reached \n', goal);
            end;
        end;
        
        % Save output values
        bestind = pops{ibest}{1};
        bestfit = fbest;
        nite = gg;
        lastpop = pops;
        lastfit = bestfits;
        
        % Stop iterating
        break;
        
    end;
    
    % Emigration
    for island=1:ni
        
        % Select origin and destination islands
        dest = island; % Destination island
        orig = island+1; if orig>ni, orig=1; end; % Origin island
        
        % Migrate individuals
        for migrant=1:nemi
            
            % Individual that will be killed by the migrant
            killed = length(pops{island}) - migrant + 1;
            
            % Migrate individual (kill the replaced individual)
            pops{dest}{killed} = pops{orig}{migrant};
            
            % Extended info
            if ninfo>10
                fprintf(['aga_islands emigrates from island=%d ',...
                    'indi=%d to replace island=%d indi=%d \n'],...
                    orig,migrant,dest,killed);     
            end;
            
        end;
        
    end;

end;

end

