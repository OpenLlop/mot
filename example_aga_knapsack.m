%% Example AGA
% Find minima of a function with Genetic Algorithm (GA)
%
% Programmers:   Manel Soria         (UPC/ETSEIAT)
%                David de la Torre   (UPC/ETSEIAT)
%                Arnau Miro          (UPC/ETSEIAT)
% Date:          23/11/2016
% Revision:      4

function example_aga_knapsack
% Knapsack problem solved with GA
% In this example the search space is a set of [0,1] integer values.
%
% Each item can be inside the sack (1) or outside (0). If it is inside,
% contributes to value but also adds weight. We look for maximum value, for
% total weight <= capacity.
%
% An individual will be a vector of 0,1 values, represented using floating
% point numbers for simplicity, but aga would allow different data types.

% Define heuristic function options (optional)
opts.ninfo = 2; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 1; % Parallel execution of fitness function
opts.nhist = 1; % Saved history level (0=none, 1=fitness, 2=full)

% Define Knapsack parameters
rng(50); % Use a known random seed for repetabiliy of the problem
NN = 64; % Number of items in the knapsack
maxval = 100; % maximum possible value of one item
maxw = 20; % max weight 
weights = randi(maxw,NN,1); % Weights
values = randi(maxval,NN,1); % Values
capacity = floor(sum(weights)/2); % Capacity
HUGE = NN*maxval*100; % Huge number

% Define AGA parameters
goal = -Inf; % Target fitness value
ng = 100; % Number of generations
np = 100; % Population size
N = [5,... % Number of elites
    floor(np*0.4),... % Number of mutants
    floor(np*0.05),...% Number of newcomers
    floor(np*0.1)]; % Number of parents

% Randomize random seed
rng('shuffle'); % We don't want repetability in the heuristic 

% Now we construct an initial population list
pop = cell(1,np);
for i=1:np
    pop{i} = ranfun(); 
end

% Execute Genetic Algorithm
[bestind, bestfit, nite, lastpop, lastfit, history] = aga ( ...
    opts, pop, goal, ng, N, @unifun, @fitfun, ...
    @mutfun, @repfun, @ranfun, @prifun );

                    
%% Plot history

% Number of generations in history array
ngens = size(history,1);

% Get fitness history
if opts.nhist>1 && iscell(history) % Full history; get fitness values
    fithist = zeros(ngens,1);
    for i=1:ngens
        fithist(i) = history{i,2}(1);
    end
else % Simple history
    fithist = history;
end

% Create figure
figure('Position',[400,200,900,600]);

% Plot fitness history
plot(fithist,'o-');

% Beautify plot
grid minor;
title('Genetic Algorithm optimization | Knapsack problem');
xlabel('Generation [#]');
ylabel('Best fitness function value');


%% Functions

    % Random individual
    function x = ranfun()
        x = transpose(randi([0,1],[1,NN]));
    end

    % Fitness function
    function f = fitfun(x) 
        f = -dot(x,values)+HUGE*NN*(dot(x,weights)>capacity);
    end

    % Discard identical individuals
    function [pop,f] = unifun(pop,f)
        % Convert population to a matrix, individuals by rows
        %mat = transpose(reshape(cell2mat(pop),NN,np));
        %matu = transpose(unique(mat,'rows')); % Get unique individuals
        %popu = mat2cell(matu,NN,ones(1,size(matu,2)));
        %f = unique(f);
    end

    % Mutate an individual
    function x = mutfun(x,f) 
        for ii=1:4
            p=randi(size(x,1));
            x(p)=1-x(p);
        end
    end
                
    % Reproduction of two individuals
    function c = repfun(a,b,fa,fb)
        e=size(a,1);
        p=randi(size(a,1));
        c=[a(1:p-1);b(p:e)];
    end

    % Print and individual
    function prifun(x)
        for ii=1:size(x,1)
            fprintf('%d',x(ii,1));
        end
    end

end
