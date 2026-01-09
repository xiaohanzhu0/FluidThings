function save_all(x1, x2, res_list, cost_list)
%------------------------------------------------------------
% 1) make a folder named by the current date/time
%------------------------------------------------------------
timestamp   = datestr(now, 'yyyymmdd_HHMMSS');   % e.g. '20250623_143015'
folderName  = ['./outputs/' timestamp];          % prepend any label you like
mkdir(folderName);                               % create it

%------------------------------------------------------------
% 2) find all open figures and save them
%------------------------------------------------------------
figs = findall(0, 'Type', 'figure');
for k = 1:numel(figs)
    fig = figs(k);
    
    % build filenames
    pngFile = fullfile(folderName, sprintf('Figure_%02d.png', k));
    figFile = fullfile(folderName, sprintf('Figure_%02d.fig', k));  % optional
    
    % save as PNG
    saveas(fig, pngFile);
    
    % also save as MATLAB .fig (uncomment if you want)
    % savefig(fig, figFile);
end

%------------------------------------------------------------
% 3) copy your script into the new folder
%------------------------------------------------------------
srcScript = mfilename('fullpath');   % full path of this running file
[~, name, ext] = fileparts(srcScript);
if strcmp(name, 'myrun')  % if this is myrun.m
    copyfile(srcScript, fullfile(folderName, ['myrun' ext]));
else
    % if running from another file, adjust to copy myrun.m explicitly:
    copyfile('myrun.m', fullfile(folderName, 'myrun.m'));
end

save(fullfile(folderName, 'result.mat'),"x1","x2","res_list", "cost_list");

end