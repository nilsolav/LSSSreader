%% clean slate
clear
close all
clc

%% Get the paths to the test data

whr = which('LSSSreader_readsnapfiles');
[dr,~,~] = fileparts(whr);

par(1).raw_file  = fullfile(dr,'/../exampledata/S2006101/tokt2006101-D20060124-T030844.raw');
par(1).snap_file = fullfile(dr,'/../exampledata/S2006101/tokt2006101-D20060124-T030844.work');

par(2).raw_file  = fullfile(dr,'/../exampledata/S2016837/2016837-D20160503-T093515.raw');
par(2).snap_file = fullfile(dr,'/../exampledata/S2016837/2016837-D20160503-T093515.work');

par(3).raw_file  = fullfile(dr,'/../exampledata/S2005114/tokt2005114-D20051118-T062010.raw');
par(3).snap_file  = fullfile(dr,'/../exampledata/S2005114/tokt2005114-D20051118-T062010.snap');


for i=1:length(par)
    if ~exist(par(i).raw_file)|~exist(par(i).snap_file)
        disp('Files are not present')
    end
end

%% Pick a file
file=3;
snap = par(file).snap_file;
raw  = par(file).raw_file;

%% Read snap file
layer = LSSSreader_readsnapfiles(snap);

%% Read raw file and convert to sv
[raw_header,raw_data] = readEKRaw(raw);
raw_cal = readEKRaw_GetCalParms(raw_header, raw_data);
Sv = readEKRaw_Power2Sv(raw_data,raw_cal);

% Get the transducer depth
f=1;
td = double(median(raw_data.pings(f).transducerdepth));

%% Plot result
[fh, ih] = readEKRaw_SimpleEchogram(Sv.pings(f).Sv, 1:length(Sv.pings(f).time), Sv.pings(f).range);

%% Plot the interpretation mask
hold on
cs = pink;

for i=1:length(layer)
    %plot(layer(i).x,layer(i).y)
    patch(layer(i).x,layer(i).y-td,cs(i,:),'FaceColor',cs(i,:),'FaceAlpha',.1)
    
end




