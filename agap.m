function [ bestind, bestfit, nite, lastpop, lastfit, history ] = agap ( ...
    opts, pop, goal, ng, N, npf, unifun, fitfun, ...
    mutfun, repfun, ranfun, prifun )
%AGA finds minimum of a function using Genetic Algorithm (GA)
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre   (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          10/05/2018
%Revision:      4
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
%   goal:       if function value is below goal, iterations are stopped
%   ng:         maximum number of generations allowed
%   N:          population control parameters
%       N(1)        ne: number of elite individuals that remain unchanged
%       N(2)        nm: number of mutants
%       N(3)        nn: number of newcomers
%                   The rest are descendants
%       N(4)        na: number of parents. The descendants are choosen
%                   among the na best individuals
%                   The rest of individuals (ie: nn=length(pop)-ne+nm+nd)
%                   are newcomers, randomly choosen
% 
%   If there are less than nm-1 non-identical indivials, population is 
%   considered degenerate and iterations stop
%
%   Call back functions to be provided by the user:
%   unifun:     deletes repeated individuals in a population
%               receives a population+fitness and returns a
%               population+fitness (a population is a list of individuals)
%   fitfun:     fitness function, given one individual returns its fitness
%               (RECALL that in this GA algoritm fitness is MINIMIZED)
%               fitness vector is: [fit_obj1, fit_obj2, ..., fit_objN]
%   mutfun:     mutation funcion, given one individual and its fitness,
%               mutfun should return a mutant individual. Fitness is given
%               in case mutation intensity is to be decreased when close
%               to the goal
%   repfun:     given two individuals and their fitnesses, returns a
%               descendant
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

% Get configuration options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else, ninfo = 1; end
if isfield(opts,'label'), label = opts.label; else, label = 0; end
if isfield(opts,'dopar'), dopar = opts.dopar; else, dopar = 0; end
if isfield(opts,'nhist'), nhist = opts.nhist; else, nhist = 1; end

% Create output folder for simulation frames
date_now = datetime('now','Format','yyyyMMddHHmmss');
dirname = fullfile('output','pareto',string(date_now));
mkdir(dirname);

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
nn = N(3); % Number of newcomers
na = N(4); % Number of breeders (selected from the best individuals)
np = length(pop); % Population size
nd = np - N(1) - N(2) - N(3); % Number of descendants

% Check number of objectives of fitfun
fi0 = feval(fitfun,pop{1});
no = size(fi0(:), 1);

% Safety checks
if na<=0, na=1; end
if nn<0, error('AGA nn (number of newcomers) must be positive'); end

% Iterate through generations
for g=1:ng
    
    % Preallocate variables
    fi = zeros(np, no); % Preallocate/clear fitness array

    % Clean population: remove repeated individuals
    [pop,fi] = unifun(pop,fi); % Return unique individuals
    pop = pop(~cellfun('isempty',pop)); % Remove empty individuals
    fi = fi(~cellfun('isempty',pop),:); % Remove empty fitness
    ncleanpop = length(pop); % Length of clean population
       
    % Avoid population degeneration (i.e., poor genetic pool)
    if ncleanpop<na % Clean population size is less than breeders size
        
        % Show info
        if ninfo>0
            fprintf('AGA label=%d degenerate population, leaving\n',label);
        end
        
        % Save last iteration data
        bestind = pop{1}; % Save best individual
        bestfit = fi(1,:); % Save fitness level of last best individual
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
        parfor i=1:np, fi(i,:) = feval(fitfun,pop{i}); end
    else % Serial execution
        for i=1:np, fi(i,:) = feval(fitfun,pop{i}); end
    end

    % Sort population into pareto fronts
    idx_front = sort_pareto();

    % Plot population & pareto fronts
    plot_population();
    plot_pareto_fronts();

    % Compute crowding distance between individuals on each front
    crowding_dist = sort_crowding(fi, idx_front, npf);

    % Combine Pareto front and crowding distance for sorting
    % Sort first by pareto fronts (increasing), then:
    % For each pareto group, sort by crowding distance (decreasing)
    [~, idx] = sortrows([idx_front, -crowding_dist]);

    % Sort population individuals by pareto front & crowding distance
    fi = fi(idx,:); % Sort fitness
    pop = pop(idx); % Sort population

    % Save history
    if nhist>1 % Save full history {population,fitness}
        history{g,1} = pop; %#ok
        history{g,2} = fi; %#ok
        history{g,3} = idx_front; %#ok
    elseif nhist>0 % Save best fitness only
        history{g} = fi(idx_front==1,:); %#ok
    end
    
    % Check if reached target fitness or max generations 
    if any(fi(1,:)<=goal) || g>=ng % Target achieved
        
        % Save last iteration data
        bestind = pop{1}; % Save best individual
        bestfit = fi(1,:); % Save fitness level of last best individual
        nite = g; % Save current generation index
        lastpop = pop; % Save last population
        lastfit = fi; % Save last fitness values
        
        % Show info
        if ninfo>0
            fprintf('AGA label=%d nite=%2d fitbest=%f',label,nite,bestfit);
            if ~isempty(prifun), fprintf(' best='); prifun(bestind); end
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
        fprintf('AGA label=%d g=%2d fitbest=%f',label,g,fi(1));
        if ~isempty(prifun), fprintf(' best='); prifun(pop{1}); end
        fprintf('\n');
    end
    
    % Preallocate population for next generation
    % <<[elites, mutants, descendants, newcomers]<<
    nextpop = cell(1, ne + nm + nd + nn);

    % Elites selection: top 'ne' individuals
    for i=1:ne
        nextpop{i} = pop{i}; % Copy elite into next generation
    end

    % Mutation: select 'nm' individuals based on crowding distance
    for i=1:nm
        k = ne + i; % Next individual in the population
        if isempty(mutfun), nextpop{k} = pop{k}; % Do not mutate
        else, nextpop{k} = mutfun(pop{k}, fi(k)); % Mutate
        end
    end

    % Descendants: perform crossover among selected parents
    % Parents are selected based on their ranking and crowding distance
    for i=1:nd
        k = ne + nm + i; % Next individual in the population
        parentA = randi([1,na]); % Parent A is choosen among the na best 
        parentB = randi([1,na]); % Parent B is choosen among the na best 
        nextpop{k} = repfun(pop{parentA}, pop{parentB}, ...
            fi(parentA), fi(parentB)); % Breed individuals A and B
    end

    % Newcomers: randomly generated individuals
    for i=1:nn
        k = ne + nm + nd + i; % Next individual in the population
        nextpop{k} = ranfun(); % Generate random individual
    end
    
    % Update population for the next generation
    pop = nextpop;
    
end


    % Compute dominance on all population
    function [idx_front] = sort_pareto()
    
        % Pareto fronts
        idx_front = NaN(np,1);

        % Build pareto fronts
        for pf=1:npf % for each pareto front...
    
            % Info
            if np > 1E4, fprintf("Pareto %d...\n", pf); end
    
            % Local pareto front copy (for parfor loop)
            idx_front_current = idx_front;

            % Run in parallel
            parfor ii=1:np % For each individual of population...
    
                % If individual is already assigned to front, skip
                if ~isnan(idx_front_current(ii))
                    continue
                end
            
                % Check dominance on current individual
                fi_pp = fi(ii,:); % Fitness of current individual
                dominated_pp = ones(np,1); % Assume pp is dominated
                for oo=1:no % For each objective...
            
                    % Check if pp is dominated by someone, on objective oo
                    % someone is dominant over pp: it has lower fitness
                    dominated_oo = fi(:,oo) < fi_pp(oo); %#ok
            
                    % Ignore individuals on other pareto fronts 
                    % (but not those on the current one)
                    dominated_oo(idx_front_current<pf) = 0; %#ok
            
                    % Aggregate dominance over all objectives
                    dominated_pp = dominated_pp & dominated_oo; 
            
                end
            
                % If pp is not dominated by anyone else --> pp is dominant
                pp_is_dominant = ~any(dominated_pp);
            
                % Debug info
                % fprintf("  %d: dominant? %d\n", i, pp_is_dominant);
    
                % If pp is dominant, add it to Pareto front pf
                if pp_is_dominant
                    idx_front(ii) = pf;
                end
    
            end
    
            % If all individuals are assigned to paretos, stop searching
            if ~any(isnan(idx_front))
                break;
            end
    
        end
    
    end

    % Compute crowding distance on all population on each pareto front
    function [crowding_dist] = sort_crowding(fi, idx_front, npf)

        % Initialize crowding distance array
        crowding_dist = zeros(size(fi, 1), 1);
    
        % For each pareto front...
        for pf = 1:npf

            % Find the individuals in the current Pareto front
            idy = find(idx_front == pf);
            if isempty(idy) % Skip if no individuals in this front
                continue;
            end
            
            % Extract fitness values of individuals in the current front
            fi_pf = fi(idy, :);
            
            % Initialize crowding distance for individuals in this front
            cd_pf = zeros(length(idy), 1);
            
            % Number of objectives
            n_obj = size(fi, 2);
            
            % For each objective...
            for m = 1:n_obj

                % Sort individuals based on the m-th objective
                [fi_s, idk] = sort(fi_pf(:, m));
                
                % Assign infinite crowding distance to boundary individuals
                cd_pf(idk(1)) = Inf;
                cd_pf(idk(end)) = Inf;
                
                % Maximum and minimum fitness values for normalization
                fi_max = max(fi_s);
                fi_min = min(fi_s);
                if fi_max == fi_min
                    fi_max = fi_min + 1; % Avoid division by zero
                end
                
                % Calculate crowding distances for intermediate individuals
                for ii = 2:(length(idk)-1)
                    cd_ii = (fi_s(ii+1) - fi_s(ii-1)) / (fi_max - fi_min);
                    cd_pf(idk(ii)) = cd_pf(idk(ii)) + cd_ii;
                end
            end
            
            % Assign calculated distances back to the main array
            crowding_dist(idy) = cd_pf;

        end

    end

    % Plot pareto fronts
    function plot_pareto_fronts()
    
        % Create figure
        fh = figure();
        fh.Position = [400,200,900,600];
        hold on;
        grid on;
        box on;
        axis equal;
        view(40,0);

        % Limits
        xlim([0,5]);
        ylim([0,5]);
        zlim([0,5]);

        % Labels
        xlabel('Objective #1: dist to xa');
        ylabel('Objective #2: dist to xb');
        zlabel('Objective #3: dist to xc');
        title('Optimising distances to 3 points');

        % Colormap
        cm = jet(npf);
        colormap(cm);

        % Marker size
        if np <= 1000, mkrsz = 30; % Small population
        elseif np <= 10000, mkrsz = 15; % Medium population
        else, mkrsz = 5; % Large population
        end

        % Plot fitnesses, color by pareto fronts
        scatter3(fi(:,1),fi(:,2),fi(:,3),mkrsz,idx_front,'o','filled');

        % Plot pareto fronts lines
        if np < 100 % Low population
            for pf=1:10
                fi_idx = fi(idx_front==pf,:);
                [~, idx_sorted] = sort(fi_idx(:,1));
                fi_idx = fi_idx(idx_sorted,:);
                plot3(fi_idx(:,1),fi_idx(:,2),fi_idx(:,3), ...
                    '-','color',cm(pf,:));
            end
        end

        % Colorbar
        hcb = colorbar;
        title(hcb, '# pareto front');

        % Save figure
        fname = fullfile(dirname,sprintf('pf_gen_%03d',g));
        print(fh,'-dpng','-r300',fname + '.png');
        savefig(fh,fname + '.fig');
        close(fh);

    end

    % Plot population
    function plot_population()
    
        % Create figure
        fh = figure();
        fh.Position = [400,200,900,600];
        hold on;
        grid on;
        box on;
        axis equal;
        view(0,90);
        
        % Limits
        xlim([0,10]);
        ylim([0,10]);
        zlim([0,npf]);

        % Scale z-axis down
        daspect([1,1,npf/5]);
        
        % Labels
        xlabel('x');
        ylabel('y');
        zlabel('# pareto front');
    
        % Colormap
        cm = jet(npf);
        colormap(cm);

        % Equilateral triangle
        x0 = [5,5]'; % Offset
        xa = x0 + [0,0]';
        xb = x0 + [2,0]';
        xc = x0 + [1,sqrt(3)]';
    
        % Plot targets
        plot3(xa(1),xa(2),0,'ko');
        plot3(xb(1),xb(2),0,'ko');
        plot3(xc(1),xc(2),0,'ko');
        plot3([xa(1),xb(1)],[xa(2),xb(2)],[0,0],'k-');
        plot3([xb(1),xc(1)],[xb(2),xc(2)],[0,0],'k-');
        plot3([xa(1),xc(1)],[xa(2),xc(2)],[0,0],'k-');
        
        % Marker size
        if np <= 200, mkrsz = 5; % Smol population
        elseif np <= 1000, mkrsz = 10; % Small population
        elseif np <= 10000, mkrsz = 5; % Medium population
        else, mkrsz = 3; % Large population
        end

        % Plot individuals
        for ii=1:np
    
            % Select plotting marker
            if np < 200 % Smol population
                if ii<=ne, marker = 'v'; % Elites (rv)
                elseif ii<=ne+nm, marker = 'o'; % Mutants (mo)
                elseif ii<=ne+nm+nd, marker = 'x'; % Descendants (bx)
                else, marker = 's'; % Newcomers (ks)
                end
            else % High population
                marker = '.';
            end
            
            % Select color
            if isnan(idx_front(ii)), c = 'k';
            else, c = cm(idx_front(ii),:); end
            if np > 1000 % High population
                if idx_front(ii) == 1, c = 'r'; end
            end
    
            % Plot individual
            x = pop{ii}(1);
            y = pop{ii}(2);
            z = idx_front(ii);
            ph = plot3(x,y,z,marker,'MarkerSize',mkrsz,'color',c);
    
            % Save legend ticks
            if ii==ne, lh(1) = ph; % Elite
            elseif ii==ne+nm, lh(2) = ph; % Mutant
            elseif ii==ne+nm+nd, lh(3) = ph; % Descendant
            elseif ii==ne+nm+nd+1, lh(4) = ph; % Newcomer
            end

        end

        % Legend (only for smol population)
        if np < 200
            legend(lh(1:4),'Elites','Mutants','Descendants','Newcomers',...
                'Location','NorthEastOutside');
        end

        % Colorbar (with ticks manually adjusted)
        hcb = colorbar;
        title(hcb, '# pareto front');
        hcb.Ticks = linspace(0, 1, npf);
        hcb.TickLabels = cellstr(num2str((1:npf)'));

        % Save figure
        fname = fullfile(dirname,sprintf('pop_gen_%03d',g));
        print(fh,'-dpng','-r300',fname + '.png');
        savefig(fh,fname + '.fig');
        close(fh);

    end

end

