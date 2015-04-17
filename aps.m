function [ bestind, bestfit, nite, lastpop, lastfit, history ] = ...
    aps( opts, pop, v, c1, c2, P1, nitemax, goal, ... 
    fitfun, prifun )
%APS finds minimum of a function using Particle Swarm (PS)
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre   (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          14/04/2015
%Revision:      1
%
%Usage:         [bestind, bestfit, nite, lastpop, lastfit, history] = ...
%                   APS( opts, pop, v, c1, c2, P1, nite, goal, ...
%                   fitfun, prifun )
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
%   pop:        initial population
%   v:          initial population velocity
%   c1:         local learning factor
%   c2:         global learning factor
%   P1:         maximum step size (velocity) allowed in one iteration
%               expressed as a fraction of the domain size
%               The domain size is estimated based on the initial
%               distribution of the particles
%               P1 should be about 1E-2
%               Low values of P1 slow down the algorithm, too high values
%               can be unstable
%   nitemax:    maximum number of iterations allowed
%   goal:       If function value is below goal, iterations are stopped
%
%   Call back functions to be provided by the user:
%   fitfun:     fitness function
%   prifun:     Prints individual
%
%Outputs:
%   bestind:    best individual from the last generation
%   bestfit:    fitness value of best individual from the last population
%   nite:       number of iterations performed
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

% Preprocessing
np = length(pop); % Population size
fi = zeros(np,1); % Fitness of each particle
fib = zeros(np,1); % Personal best fitness of each particle
popb = pop; % Personal best position of each particle
bestfit = 0; % Best global fitness

% Estimate the domain size to limit particle velocity
xmax = pop{1}; % Select first particle as starter
xmin = pop{1};
for i=2:np % Get the global max/min position of the entire population
    xmax = max(xmax,pop{i}); % Global maximum position
    xmin = min(xmin,pop{i}); % Global minimum position
end;
vmax = P1 * norm(xmax-xmin); % Maximum allowed velocity

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
    
    % Get global minimum fitness value and its particle
    [minfit, mfp] = min(fi);

    % If a better particle is found, update global best values
    if minfit<bestfit || ite==1
        bestind = pop{mfp}; % New global best particle
        bestfit = minfit; % New global best fitness
    end;
    
    % Save history
    if nhist>1 % Save full history {population,fitness}
        history{ite,1} = pop; %#ok
        history{ite,2} = fi; %#ok
    elseif nhist>0 % Save best fitness only
        history(ite) = fi(1); %#ok
    end;
    
    % Show info if required
    if ninfo>1
        fprintf('PS label=%d ite=%d',label,ite); % Print info
        fprintf(' best= '); prifun(pop{mfp}); % Print best individual
        fprintf(' vbest= '); prifun(v{mfp}); % Print best individual v
        fprintf(' fbest=%e \n',minfit); % Print best individual fitness
    end;
    
    % Check if reached target fitness or max iterations
    if minfit<goal || ite>=nitemax % Target achieved
        
        % Save last iteration data
        bestind = pop{mfp}; % Save best individual
        bestfit = minfit; % Save fitness level of last best individual
        nite = ite; % Save current generation index
        lastpop = pop; % Save last population
        lastfit = fi; % Save last fitness values
        
        % Stop iterating
        break;
        
    end;
    
    % Update positions
    for i=1:np
        pop{i} = pop{i} + v{i}; % Add velocity for one time step
    end

    % Update velocities
    for i=1:np
        v{i} = v{i} ... % Previous velocity
            + c1 * rand() * (popb{i} - pop{i}) ... % Local learning 
            + c2 * rand() * (bestind - pop{i}); % Global learning
    end;
   
    % Limit velocities
    for i=1:np
        nvi = norm(v{i}); % Velocity norm
        if nvi>vmax % Ecessive velocity
            v{i} = v{i} * vmax / nvi; % Scale velocity
        end;
    end;
    
end;

end

