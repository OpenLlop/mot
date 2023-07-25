function [ ] = savefun_rastrigin ( pop, fi, g )

% Get the current process PID
try % For Windows
    pid = feature('GetPid');
catch
    try % For UNIX
        [~, result] = system('echo $$');
        pid = str2double(result);
    catch
        pid = "UNK";
    end
end

% Backup dir & filename
odir = "backups/";
casename = 'rastrigin';
fname = sprintf('%s_%d_%d', casename, g, pid);

% Create backup dir if required
if ~exist(odir, 'dir'), mkdir(odir); end

% Save vars into mat file
save(odir + fname + ".mat", "pop", "fi", "g");

end

