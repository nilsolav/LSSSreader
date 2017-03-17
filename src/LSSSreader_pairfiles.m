function directory=LSSSreader_pairfiles(directory)
%
% this function lists pairs of snap, work and raw files from multiple
% directories. If only the snapdir is given it assumes that all files are
% in the same directory.
%
% Input:
% directory(i)          : Directory j where the files are located
% directory(i).snap_dir : The location of the snap files
% directory(i).work_dir : The location of the work files (optional)
% directory(i).raw_dir  : The location of the raw files (optional)
%
% Output:
% directory(i).files{j,1}  : full path to snapfile
% directory(i).files{j,2}  : full path to workfil
% directory(i).files{j,3}  : full path EK raw file

% Loop over direcotries
for i=1:length(directory)
    
    % Get the snapfile
    snap = dir(fullfile(directory(i).snap_dir,'*.snap'));
    
    % Get the workfiles
    if isfield('work_dir',directory)
        work = dir(fullfile(directory(i).work_dir,'*.work'));
    else
        work = dir(fullfile(directory(i).snap_dir,'*.work'));
    end
    
    % Get the raw files
    if isfield('raw_dir',directory)
        raw = dir(fullfile(directory(i).raw_dir,'*.raw'));
    else
        raw = dir(fullfile(directory(i).snap_dir,'*.raw'));
    end
    
    % Combine the files
    k = 1;
    N = max([length(raw) length(snap) length(work)]);
    files={};
    files=combinefiles(files,snap,1);
    files=combinefiles(files,work,2);
    files=combinefiles(files,raw,3);
    directory(i).files = files;
end


function files0 = combinefiles(files0,files,fp)
% List
s=size(files0);
for i=1:length(files)
    file = fullfile(files(i).folder,files(i).name);
    [~,filemain,~] = fileparts(file);
    % does the file already exist?
    ex = false;
    for j=1:s(1)
        for l=1:s(2)
            if ~isempty(files0{j,l})
            [~,filemain0,~] = fileparts(files0{j,l});
            if strcmp(filemain,filemain0)
                ex=true;
                files0{j,fp} =file;
            end
            end
        end
    end
    if ~ex
        % Append the file
        files0{s(1)+1,fp} = file;
        s(1)=s(1)+1;
    end
end


