%% Initialize

whr = which('LSSSreader_readsnapfiles');
[dr,~,~] = fileparts(whr);
dr = dr(1:end-3);

% Location of example files for Gavin
%dr = 'F:\Documents\Dropbox\toRoland';

% Recursively list relevant data in the example directory. This can be used
% to search files in callisto.
files.snap = rdir(fullfile(dr,'exampledata','**','*.snap'));
files.work = rdir(fullfile(dr,'exampledata','**','*.work'));
files.raw = rdir(fullfile(dr,'exampledata','**','*.raw'));

% Match the different files
files=LSSSreader_pairfiles(files);


%% Pick a file
for file=2%[1:4 6:size(files.F,1)] % I think we should remove the macseis data from the test data. It is way to big.
    snap = files.F{file,1};
    work = files.F{file,2};
    raw  = files.F{file,3};
    if isempty(snap)
        snap=work;
    end
    % Read snap file
    [school,layer,exclude,erased] = LSSSreader_readsnapfiles(snap);
    if true    % Set to false if you do not wan't to plot the echograms
        % Read raw file and convert to sv
        [raw_header,raw_data] = readEKRaw(raw);
        raw_cal = readEKRaw_GetCalParms(raw_header, raw_data);
        Sv = readEKRaw_Power2Sv(raw_data,raw_cal);
        
        % Get the transducer depth
        if length(raw_data.pings) > 1
            ch=2; % Use 38 kHz if we can (which is usually channel 2 on IMR ships)
            f='38';
        else
            ch=1;
            f='38';
        end
        td = double(median(raw_data.pings(ch).transducerdepth));
        
        % Plot result
        for ch =1:length(Sv.pings)
            [fh, ih] = readEKRaw_SimpleEchogram(Sv.pings(ch).Sv, 1:length(Sv.pings(ch).time), Sv.pings(ch).range);
            hold on
            % Plot the interpretation layer
            LSSSreader_plotsnapfiles(layer,school,erased,exclude,num2str(Sv.pings(ch).frequency(1)/1000),td)
        end
    else
        f='38';
        td=0;
        ch=2;
        figure
        maxRange = [];
        hold on
        % Plot the interpretation layer
        LSSSreader_plotsnapfiles(layer,school,erased,exclude,f,td)
    end
end
