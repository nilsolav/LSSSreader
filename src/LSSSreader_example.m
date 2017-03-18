%% clean slate for some of us
user = char(java.lang.System.getProperty('user.name'));
if ~strcmp(user, 'gavinj')
    clear
    close all
    clc
else
    warning('I must ensure that I do not amuse my self too much. Not healthy.')
end

whr = which('LSSSreader_readsnapfiles');
[dr,~,~] = fileparts(whr);

directory(1).snap_dir = fullfile(dr,'/../exampledata/S2006101/');
directory(2).snap_dir = fullfile(dr,'/../exampledata/S2016837/');
directory(3).snap_dir = fullfile(dr,'/../exampledata/S2005114/');

d=LSSSreader_pairfiles(directory);


%% Get the paths to the test data


par(1).raw_file  = fullfile(dr,'/../exampledata/S2006101/tokt2006101-D20060124-T030844.raw');
par(1).snap_file = fullfile(dr,'/../exampledata/S2006101/tokt2006101-D20060124-T030844.work');

par(2).raw_file  = fullfile(dr,'/../exampledata/S2016837/2016837-D20160503-T093515.raw');
par(2).snap_file = fullfile(dr,'/../exampledata/S2016837/2016837-D20160503-T093515.work');

par(3).raw_file  = fullfile(dr,'/../exampledata/S2005114/tokt2005114-D20051118-T062010.raw');
par(3).snap_file  = fullfile(dr,'/../exampledata/S2005114/tokt2005114-D20051118-T062010.snap');

% Example with two types of erased region and an exclude region
par(4).raw_file  = fullfile(dr,'/../exampledata/S2014119/2014119-D20141029-T172311.raw');
par(4).snap_file  = fullfile(dr,'/../exampledata/S2014119/2014119-D20141029-T172311.work');


for i=1:length(par)
    if ~exist(par(i).raw_file)|~exist(par(i).snap_file)
        disp('Files are not present')
    end
end

%% Pick a file
for file=1:4
snap = par(file).snap_file;
raw  = par(file).raw_file;

%% Read snap file
[school,layer,exclude,erased] = LSSSreader_readsnapfiles(snap);

%% Read raw file and convert to sv
[raw_header,raw_data] = readEKRaw(raw);
raw_cal = readEKRaw_GetCalParms(raw_header, raw_data);
Sv = readEKRaw_Power2Sv(raw_data,raw_cal);

% Get the transducer depth
f=1;
if length(raw_data.pings) > 1
    f=2; % Use 38 kHz if we can (which is usually channel 2 on IMR ships)
end
td = double(median(raw_data.pings(f).transducerdepth));

%% Plot result
[fh, ih] = readEKRaw_SimpleEchogram(Sv.pings(f).Sv, 1:length(Sv.pings(f).time), Sv.pings(f).range);

% Plot the interpretation mask
hold on
cs = cool;
for i=1:length(layer)
    if length(layer)>1
        col = round(interp1(linspace(1,length(layer),size(cs,1)),1:size(cs,1),i));
    else
        col=1;
    end
    patch(layer(i).x,layer(i).y-td,cs(col,:),'FaceColor',cs(col,:),'FaceAlpha',.3)
end

cs = hot;
for i=1:length(school)
    col = round(interp1(linspace(1,length(school),size(cs,1)),1:size(cs,1),i));
    patch(school(i).x,school(i).y-td,cs(col,:),'FaceColor',cs(col,:),'FaceAlpha',.3)
end


end
