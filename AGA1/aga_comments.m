function [lop,fop,nite]=aga_comments(ninfo,label, ... 
                           pop, ... 
                           ng,nm,nr,nn, goal, ... 
                           funique,fitfun,mutfun,reprofun,ranfun,prifun)  
% Iterates to find mimumum of a function using Genetic Algorithm
% (c) 2013 - Manel Soria - ETSEIAT
%
% ninfo:    iteration control; prints every ninfo iterations 
% label:    integer number that precedes the prints in case output is to be
%           filtered
% pop:      list with initial population elements
% ng:       number of generations
% nm:       control of elite individuals: after sorting by fitness, 
%           pop{1} .. pop{nm-1} remain unchanged
% nr:       control of mutations: pop{nm} to pop{nr-1} become mutations 
%           of the elite 
% nn:       control of reproduction: pop{nr}..pop{nn-1} are descendants of 
%           random individuals selected from all the population
%           the remaining pop{nn}..pop{end} are newcommers, random 
%           individuals 
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
% reprofun: Given two individuals, returns a descendant 
% ranfun:   Returns a random individual
% prifun:   Prints individual
%
% aga returns:
% lop:      list with the population sorted by  fitness  
% fop:      minimum value of fitfun found
% nite:     number of iterations performed 
               
% Size of initial population (length of initial population array)
np = length(pop);

% Iterate through generations
for g=1:ng
    
    % Save current generation index
    nite = g;
    
    % Clean population: delete repeated individuals
    % This is done to avoid useless fitness func evaluations
    pop = funique(pop);
    
    % Avoid population degeneration (=poor genetic pool)
    if length(pop) < nm-1 % Cleaned population size less than req. elites
        fprintf('GA label=%d: degenerated population\n',label);
        break; % Break GA execution
    end;
    
    % Repopulation: fill cleaned population pool with new breedings
    for i=length(pop)+1:np % Refill population up to initial size
        pop{i} = ranfun(); % Create new "random" individual
    end;
    
    % Evaluate fitness function
    parfor i=1:np % Parallelized for perforance
        fi(i) = feval(fitfun,pop{i}); % Evaluate fitness of individual i
    end;
    
    % Sort population individuals by their fitness level
    [fi,i] = sort(fi); % Sort fitness by increasing value (lower is best)
    pop = pop(i); % Sort population by their fitness value

    % Show info
    if ninfo~=0 && mod(g,ninfo)==0
        fprintf('GA label=%d g=%3d ng=%d nm=%d nr=%d best=%e ',...
            label,g,ng,nm,nr,fi(1)); % Print general information
        if ~isempty(prifun), prifun(pop{1}); end; % Show first individual
        fprintf('\n'); % Print newline
    end;

    % Simulation end: either reached target fitness or max generations 
    if fi(1) < goal || g == ng
        
        % Save last population data
        lop = pop; % Save entire population
        fop = fi(1); % Save best fitness level
        
        % Show info
        if ninfo~=0
            if fi(1) < goal, fprintf('GA label=%d goal=%e achieved !!\n',label,goal);
            else fprintf('GA label=%d goal=%e lograt=%e maxIter, leaving\n',label,goal,fop);
            end;
        end;

        % Return from function
        return;
        
    end;
    
    % Population will now be computed as:
    % <<[elites, mutants, offsprings, newcomers]<<
    
    % Mutation: generate mutant individuals
    for i=nm:nr-1
        pop{i} = mutfun(pop{i},fi(i)); % Mutate individual
    end;
    
    % Reproduction: generate offspring individuals
    % Use random individuals from the entire population pool as parents
    for i=nr:nn-1 
        parentA = randi([1,np]);
        parentB = randi([1,np]);
        pop{i} = reprofun(pop{parentA},pop{parentB});
    end;
    
    % Newcommers: generate new individuals up to initial population size
    for i=nn:np
        pop{i}=ranfun();
    end;

end;

% You should only reach this point if your population is degenerated

% Re-evaluate fitness function with degenerated population
parfor i=1:np
    fi(i) = feval(fitfun,pop{i});
end;

% Sort degenerated population data
[fi,i] = sort(fi);
pop = pop(i); 

% Save values
lop = pop; 
fop = fi(1);

end

