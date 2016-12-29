function [ bestind, bestfit, nite, lastpop, lastfit, history ] = aps ( ...
    opts, pop, goal, DATA )
%APS finds minimum of a function using Particle Swarm (PS)
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre   (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          17/11/2016
%Revision:      2
%
%Usage:         [bestind, bestfit, nite, lastpop, lastfit, history] = ...
%                   APS( opts, pop, goal, DATA )
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
%   pop:        initial population
%   goal:       If function value is below goal, iterations are stopped
%   DATA:       structure with the specific parameters and callback
%               functions of the APS heuristic function.
%       nitemax:maximum number of iterations allowed
%       v:      initial population velocity
%       c1:     local learning factor
%       c2:     global learning factor
%       vmax:   maximum value for particle velocity
%
%       Call back functions to be provided by the user:
%       fitfun: fitness function
%       posfun: position update function
%       velfun: velocity update function
%       vscfun: velocity scaling function
%       prifun: prints best individual
%       rpsfun: returns a random position
%       rvlfun: returns a random velocity
%
%Outputs:
%   bestind:    best individual from the last generation
%   bestfit:    fitness value of best individual from the last population
%   nite:       number of iterations performed
%   lastpop:    population of last generation
%   lastfit:    fitness values of last population
%   history:    array with saved history array

% Get options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else, ninfo = 1; end;
if isfield(opts,'label'), label = opts.label; else, label = 0; end;
if isfield(opts,'dopar'), dopar = opts.dopar; else, dopar = 0; end;
if isfield(opts,'nhist'), nhist = opts.nhist; else, nhist = 1; end;

% Get heuristic parameters from data structure
nitemax = DATA.nitemax;
v = DATA.v;
c1 = DATA.c1;
c2 = DATA.c2;
vmax = DATA.vmax;
fitfun = DATA.fitfun;
posfun = DATA.posfun;
velfun = DATA.velfun;
vscfun = DATA.vscfun;
prifun = DATA.prifun;
ranfun = DATA.ranfun;
rvlfun = DATA.rvlfun;

% Build population if required
if isnumeric(pop) % Population size is given as input
    np = pop; % Population size
    pop = cell(1,np); % Preallocate population variable
    for i=1:np % Fill population
        pop{i} = ranfun(); % Generate random particle
    end;
end;

% Create history array
history = [];

% Preprocessing
np = length(pop); % Population size
fi = zeros(np,1); % Fitness of each particle
fib = zeros(np,1); % Personal best fitness of each particle
popb = pop; % Personal best position of each particle
bestfit = 0; % Best global fitness

% Build initial velocity if required
if isnumeric(v) % Population size is given as input
    vfact = v; % Get scaling factor for the velocity
    v = cell(1,np); % Preallocate velocity variable
    for i=1:np % Fill population
        v{i} = rvlfun(vfact); % Generate random velocity
        v{i} = vscfun(v{i},vmax); % Limit velocity
    end;
end;

% Iterate until convergence or max iterations
for ite=1:nitemax

    % Evaluate fitness function
    if dopar % Parallel execution
        parfor i=1:np, fi(i) = feval(fitfun,pop{i}); end;
    else % Serial execution
        for i=1:np, fi(i) = feval(fitfun,pop{i}); end;
    end;

    % Update personal best of each particle
    for i=1:np
        if fi(i)<fib(i) || ite==1
            fib(i) = fi(i); % Best fitness
            popb{i} = pop{i}; % Best particle position
        end;
    end;
    
    % Sort population individuals by their fitness level
    [fi,i] = sort(fi); % Sort fitness by increasing value (lower is best)
    pop = pop(i); % Sort population individuals
    fib = fib(i); % Sort population personal best fitness
    popb = popb(i); % Sort population personal best position
    v = v(i); % Sort population velocities
    
    % If a better particle is found, update global best values
    if fi(1)<bestfit || ite==1
        bestind = pop{1}; % New global best particle
        bestfit = fi(1); % New global best fitness
    end;
    
    % Save history
    if nhist>1 % Save full history {population,fitness}
        history{ite,1} = pop; %#ok
        history{ite,2} = fi; %#ok
    elseif nhist>0 % Save best fitness only
        history(ite) = fi(1); %#ok
    end;
    
    % Check if reached target fitness or max iterations
    if fi(1)<goal || ite>=nitemax % Target achieved
        
        % Save last iteration data
        bestind = pop{1}; % Save best individual
        bestfit = fi(1); % Save fitness level of last best individual
        nite = ite; % Save current generation index
        lastpop = pop; % Save last population
        lastfit = fi; % Save last fitness values
        
        % Show info
        if ninfo>0
            fprintf('APS label=%d nite=%2d fitbest=%f',label,nite,bestfit);
            if ~isempty(prifun), fprintf(' best='); prifun(bestind); end;
            if bestfit<goal % Goal achieved
                fprintf(' goal=%e achieved, leaving\n',goal);
            else % Maximum generations reached (goal not achieved)
                fprintf(' max. iterations reached, leaving\n');
            end;
        end;
        
        % Stop iterating
        break;
        
    end;
    
    % Show info if required
    if ninfo>1
        fprintf('APS label=%d ite=%2d fitbest=%e',label,ite,fi(1));
        if ~isempty(prifun), fprintf(' best='); prifun(pop{1}); end;
        fprintf('\n');
    end;
    
    % Update positions
    for i=1:np
        pop{i} = posfun(pop{i},v{i}); % Add velocity for one time step
    end

    % Update velocities
    for i=1:np
        v{i} = velfun(v{i},pop{i},popb{i},bestind,c1,c2);
    end;
   
    % Limit velocities
    for i=1:np
        v{i} = vscfun(v{i},vmax);
    end;
    
end;

end

