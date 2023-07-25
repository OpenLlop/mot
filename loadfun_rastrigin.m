function [ pop ] = loadfun_rastrigin ( g, nn )

% Specific configuration for loading backupped inds
at_very_gens = 5;

% Declare vars
pop = cell(1, nn);
ki = 1;

% Load backupped only every XX generations
if mod(g, at_very_gens) == 0

    % Info
    disp("Backupped festival!");

    % Get the current process PID
    try 
        pid = feature('GetPid'); % For Windows
    catch
        [~, result] = system('echo $$');
        pid = str2double(result); % For UNIX
    end
    
    % Output filename
    odir = "backups/";
    
    % Get the list of files and folders in the specified directory
    odir_contents = dir(odir);
    
    % Filter out the directories (if you only want files)
    file_list = {odir_contents(~[odir_contents.isdir]).name};
    
    % Defult file to load: none
    target_file_to_load = [];
    
    % Select a random file to load individuals from
    for ii=1:100
        
        % Try one random file from the backup dir 
        idx = randi(length(file_list));
        target_file_to_load = file_list{idx};
    
        % Chop-chop name
        parts = strsplit(target_file_to_load, '.');
        parts = strsplit(parts{1}, '_');
        ext_g = str2double(parts{2});
        ext_pid = str2double(parts{3});
    
        % Ensure we do not have selfcest
        if ext_pid == pid
            continue;
        end
    
        % Ensure we do not breed with the ancients
        threshold = 10;
        if ext_g < g && abs(ext_g - g) > threshold
            continue;
        end
    
        % We are now ethically allowed to breed :}
        break;
    
    end
    
    % Try to get individuals from previously saved populations
    if ~isempty(target_file_to_load)
    
        try % Try to load and process the file
            
            % Load file
            raw = load(odir + target_file_to_load);
    
            for ki=1:nn % Backupped individuals
    
                % Search for the external individual that has better fi
                pop{ki} = raw.pop{ki};
    
                % If all individuals filled, return
                if (ki >= nn) 
                    return; 
                end
    
                % If no more available bakupped individuals, stop
                if (ki >= length(raw.fi)) 
                    break; 
                end
    
            end
    
        catch ex % Error while loading or processing the file
            disp(['Error occurred while loading the file: ', ex.message]);
        end

    end

end

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a and b

% Generate random individuals, up to filling the nn requested inds
for kj=ki:nn % Newcomers
    pop{kj} = ranrange(-5,5,2); % Random individual
end

end

