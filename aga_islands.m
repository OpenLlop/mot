function [ lastpops, bestfits, bestind, bestfit, history ] = ...
    aga_islands( opts, ...
    pops, ngg, nemi, ng, N, goal, ...
    funique, fitfun, mutfun, repfun, ranfun, prifun )
% Iterates to find minimum of a function using Genetic Algorithm
% (c) 2013 - Manel Soria - ETSEIAT - v1.01
% (c) 2015 - Manel Soria, David de la Torre - ETSEIAT - v1.02
%
% opts:     function control parameters [struct]
%   ninfo:  verbosity level (0=none, 1=minimal, 2=extended)
%   label:  integer number that precedes the prints in case output is to
%           be filtered
%   paral:  parallel execution of fitness function [1,0]
%   nhist:  save history (0=none, 1=fitness, 2=all{pop,fit})
%               0: history = []
%               1: history(ngg,ni) = bestfit(i)
%               2: history{ngg,ni} = [hist,nite] (refer to aga)
%   plotf:  plot fitness (0=none, 1=plot, 2=plot+save)
%               0: no plot
%               1: plot on each global generation
%               2: plot on each global generation + save plot to file
%   plotp:  plot population (0=none, 1=plot, 2=plot+save)
%               0: no plot
%               1: plot on each global generation
%               2: plot on each global generation + save plot to file
% pop:      list with initial population elements
% ngg:      number of global generations
% nemi:     number of emigrations between islands per global generation
% ng:       number of generations for each island per global generation
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
% aga_islands returns:
% lastpops: list with last populations sorted by fitness (all islands)
% lastfits: best last fitness value (all islands)
% bestind:  best individual (among all the islands)
% bestfit:  fitness value of best individual
% history:  array (ngg,ni) with the best value found after each iteration

% Set default options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else ninfo = 0; end;
if isfield(opts,'label'), label = opts.label; else label = 0; end;
if isfield(opts,'nhist'), nhist = opts.nhist; else nhist = 0; end;
if isfield(opts,'plotf'), plotf = opts.plotf; else plotf = 0; end;
if isfield(opts,'plotp'), plotp = opts.plotp; else plotp = 0; end;

% Declare history var
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
        
        % AGA options
        optsi.ninfo = ninfo-1;
        optsi.label = label + gg*1000 + island;
        optsi.nhist = nhist;
        
        % Execute AGA (Genetic Algorithm)
        [lastpop,bestfit,nite,hist] = aga( optsi, ...
            pops{island}, ng, N, goal,...
            funique, fitfun, mutfun, repfun, ranfun, prifun);
        
        % Save values at the end of local island iteration
        pops{island} = lastpop; % Save last population
        bestfits(island) = bestfit; % Save best fitness value

        % Save history
        if nhist>1 % Save full history for each island
            history{gg,island} = {hist,nite}; %#ok
        elseif nhist>0 % Save best fitness only
            history(gg,island) = bestfit; %#ok
        end;

        % Show extended info
        if ninfo>1
            fprintf(['aga_islands label= %d end ite global= %d ',...
                'illa= %d fbest= %e'],label,gg,island,bestfits(island));
            if ~isempty(prifun) % Print best individual
                fprintf(' best: '); prifun(lastpop{1});
            end;
            fprintf('\n');
        end;
        
    end;
    
    % Find best individual of all islands and island where it lives
    [fbest,ibest] = min(bestfits);
    if ninfo>0 % Show info
        fprintf('aga_islands label= %d gg= %d fbest= %8.3e island= %d ',...
            label,gg,fbest,ibest);
        if ~isempty(prifun) 
            fprintf('best: '); prifun(lastpop{i});
        end;
        fprintf('\n');         
    end;

    % Save history
    if nhist>1 % Save full history for each island
        history{gg,ni+1} = pops{ibest}{1}; %#ok
        history{gg,ni+2} = fbest; %#ok
    end;

    % Plot fitness if required
    if plotf>0
       if ~exist('fhf','var'), [fhf] = InitializePlotFitFigure(); end;
       PlotFitnessHistory(fhf); % Plot fitness history
    end;

    % Plot generation if required
    if plotp>0
       if ~exist('fhp','var'), [fhp,ph] = InitializePlotGenFigure(); end;
       PlotCurrentPopulation(fhp); % Plot generation
    end;

    % Check if reached target fitness or max generations 
    if fbest<=goal || gg>=ngg
        
        % Show info if required
        if ninfo>1
            fprintf('aga_islands stoping because ');
            if (gg==ngg) 
                fprintf('maximum number of global iterations %d has been reached \n', ng);
            else
                fprintf('goal %e has been reached \n', goal);
            end;
        end;
        
        % Save last values
        lastpops = pops{ibest};
        bestind = pops{ibest}{1};
        bestfit = fbest;
        
        % Break iteration
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


%% Auxiliar functions

    % Create fitness plot figure
    function [fh] = InitializePlotFitFigure()

        % Figure Sizing
        SS = get(0,'ScreenSize'); % Get User's Screen Size
        figW = 800; % Figure Width
        figH = 600; % Figure Height
        figW0 = SS(3)/2 - figW/2; % Initial x-Coordinate
        figH0 = SS(4)/2 - figH/2; % Initial y-Coordinate

        % Create figure
        fh = figure('Position',[figW0,figH0,figW,figH],...
            'PaperSizeMode','auto');
        
        % Plot settings
        hold on; % Hold figure
        box on; grid minor;
        title('Genetic Algorithm optimization + Islands');
        xlabel('Global generation [#]');
        ylabel('Best fitness function value');
        
    end

    % Plot fitness history
    function PlotFitnessHistory(fh)
        
        % Set current figure
        figure(fh);
        
        % Plot history
        plot(gg,fbest,'b*-');
        
        % Do events
        drawnow;
        
        % Save plot to file
        if plotf>1
            if ~isdir('output'), mkdir('output'); end; % Create dir
            print(fh,'-dpng','-r300',fullfile('output',...
                sprintf('AGAIS_FIT_%d',label))); % Print plot
        end;

    end

    % Create population plot figure
    function [fh,ph] = InitializePlotGenFigure()

        % Figure Sizing
        SS = get(0,'ScreenSize'); % Get User's Screen Size
        figW = 800; % Figure Width
        figH = 600; % Figure Height
        figW0 = SS(3)/2 - figW/2; % Initial x-Coordinate
        figH0 = SS(4)/2 - figH/2; % Initial y-Coordinate

        % Create figure
        fh = figure('Position',[figW0,figH0,figW,figH],...
            'PaperSizeMode','auto');
        
        % Plot settings
        hold on; % Hold figure
        view(0,90); box on; grid minor;

        % Create plot handles
        ph = cell(ni,np);
        
    end

    % Plot current generation
    function PlotCurrentPopulation(fh)

        % Set current figure
        figure(fh);
        
        % Title
        title({'Genetic Algorithm optimization + Islands';...
            sprintf('Generation %03d',gg)});
        
        % Plot islands
        for is=1:ni

            % Population size
            ne = N(1); % Number of elites
            nm = N(2); % Number of mutants
            ps = length(pops{is}); % Population size
            nd = ps - N(1) - N(2) - N(3); % Number of descendants

            % Delete old individuals
            if gg~=1, for ii=1:ps, delete(ph{is,ii}); end; end;

            % Plot new individuals
            for ii=1:ps

                % Select plotting marker
                if ii<=ne, marker = 'rv'; % Elites range
                elseif ii<=ne+nm, marker = 'mo'; % Mutants range
                elseif ii<=ne+nm+nd, marker = 'bx'; % Descendants range
                else marker = 'ks'; % Newcomers range
                end;

                % Plot individual
                x = pops{is}{ii}(1); % X
                y = pops{is}{ii}(2); % Y
                z = fi(ii); % Fitness
                ph{is,ii} = plot3(x,y,z,marker,'MarkerSize',4); % Plot 3D

                % Save legend ticks
                if ii==ne, lh(1) = ph{is,ii}; % Elite
                elseif ii==ne+nm, lh(2) = ph{is,ii}; % Mutant
                elseif ii==ne+nm+nd, lh(3) = ph{is,ii}; % Descendant
                elseif ii==ne+nm+nd+1, lh(4) = ph{is,ii}; % Newcomer
                end;

            end;
            
        end;

        % Legend
        legend(lh,'Elites','Mutants','Descendants','Newcomers',...
            'Location','NorthEastOutside');

        % Do events
        drawnow;
        
        % Save plot to file
        if plotp>1
            if ~isdir('output'), mkdir('output'); end; % Create dir
            print(fh,'-dpng','-r300',fullfile('output',...
                sprintf('AGAIS_POP_%d_%03d',label,g))); % Print plot
        end;

    end

end

