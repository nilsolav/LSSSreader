%% Initialize

% The example work files can be found at:
% https://github.com/nilsolav/LSSS-label-versioning
% git clone https://github.com/nilsolav/LSSS-label-versioning
%
% The .raw (and bot/idx) files needs to be copied manually from here:
% /nfs/delphi/431-Obsmet/nilsolav_scratch/LSSS-label-versioning (Linux) or
% /scratch/nilsolav/LSSS-label-versioning/ (Pallas) or
% \\delphi\felles\431-Obsmet\nilsolav_scratch\LSSS-label-versioning (windows)

dr = 'D:\DATA\LSSS-label-versioning';
if ~exist(dr,'dir') %Previous version had the example files bundled with the reader
    whr = which('LSSSreader_readsnapfiles');
    [dr,~,~] = fileparts(whr);
    dr = dr(1:end-3);
end

files.snap = rdir(fullfile(dr,'**','*.snap'));
files.work = rdir(fullfile(dr,'**','*.work'));
files.raw = rdir(fullfile(dr,'**','*.raw'));

% Match the corresponding snap, work and raw files (by file name)
files=LSSSreader_pairfiles(files);

%% Run the example data

pl = true; % Set to false for plotting the masks only (without background echograms)
%pl = false;

for file=1:size(files.F,1)
    snap = files.F{file,1};
    work = files.F{file,2};
    raw  = files.F{file,3};
    % If no snap files are present, use the work files instead
    if isempty(snap)
        snap=work;
    end
    % Read the snap (work) file
    [school,layer,exclude,erased] = LSSSreader_readsnapfiles(snap);
    if pl
        
        % Read raw file and convert to sv
        [raw_header,raw_data] = readEKRaw(raw);
        raw_cal = readEKRaw_GetCalParms(raw_header, raw_data);
        Sv = readEKRaw_Power2Sv(raw_data,raw_cal);
        
%         if length(raw_data.pings) > 1
%             ch=2; % Use 38 kHz if we can (which is usually channel 2 on IMR ships)
%             f='38';
%         else
%             ch=1;
%             f='38';
%         end
        % Get the transducer depth from the raw data
        
        % Plot the result for all available frequencies
        for ch = 1:length(Sv.pings)
            td = double(median(raw_data.pings(ch).transducerdepth));
            % Plot echogram
            [fh, ih] = readEKRaw_SimpleEchogram(Sv.pings(ch).Sv, 1:length(Sv.pings(ch).time), Sv.pings(ch).range);
            % Plot the interpretation layer
            hold on
            f=num2str(Sv.pings(ch).frequency(1)/1000);
            LSSSreader_plotsnapfiles(layer,school,erased,exclude,f,ch,td,Sv.pings(ch).time)
            title([f,' kHz'])
        end
    else
        % Plotting section when background data is not plotted
        f='38';
        td=0;
        ch=2;
        figure
        hold on
        % Plot the interpretation layer
        LSSSreader_plotsnapfiles(layer,school,erased,exclude,f,ch,td,Sv.pings(ch).time)
    end
end
