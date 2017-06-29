function LSSSreader_plotsnapfiles(layer,school,erased,exclude,f,td)
%LSSSreader_plotsnapfiles
%   LSSSreader_plotsnapfiles(layer,school,erased,exclude,ch,td) Plot layers,
%   schools, erased and exclude regions as pathces.
%
%   'layer'   -  The layer from LSSSreader_readsnapfiles
%   'school'  -  The school from LSSSreader_readsnapfiles
%   'erased'  -  The erased regions from LSSSreader_readsnapfiles
%   'exclude' -  The excluded regions from LSSSreader_readsnapfiles
%   td - The trasnducerdpeth (default 0)
%   f - the frequency to plot
%

if nargin<6
    td=0;
end

% True if the ID should be plotted
plID=true;

if ~isempty(layer)
    cs = cool;
    for i=1:length(layer)
        
        % Plot the patch
        if length(layer)>1
            col = round(interp1(linspace(1,length(layer),size(cs,1)),1:size(cs,1),i));
        else
            col=1;
        end
        patch(layer(i).x,layer(i).y-td,cs(col,:),'FaceColor',cs(col,:),'FaceAlpha',.3)

        % Plot the Id string
        if isfield(layer(i),'channel') && isempty(layer(i).channel)
            Idstring='No species ID';
        else
            Idstring=[];
            if isfield(layer(i),'channel')
                for ch = 1:length(layer(i).channel)
                    if isfield(layer(i).channel(ch),'species')
                        for sp=1:length(layer(i).channel(ch).species)
                            if strcmp(layer(i).channel(ch).frequency,f)
                                Idstring =[Idstring, ['ID:',layer(i).channel(ch).species(sp).speciesID,' fraction:',layer(i).channel(ch).species(sp).fraction,';']];
                            end
                        end
                    end
                end
            end
        end
        if isempty(Idstring)
            Idstring='No species ID';
        end
        text(layer(i).x(1),layer(i).y(1)-td,Idstring)
    end
end

Idstring=[];
% Plot the interpretation school
if ~isempty(school)
    cs = hot;
    % Loop over schools
    for i=1:length(school)
        % Plot only non empty schools (since we do not know whether an
        % empty school is assiciated to a frequency)
        if ~isempty(school(i).channel)
            % Loop over channels 
            for ch = 1:length(school(i).channel)
                % Plot only the relevant frequency
                if strcmp(school(i).channel(ch).frequency,f)
                    col = round(interp1(linspace(1,length(school),size(cs,1)),1:size(cs,1),i));
                    patch(school(i).x,school(i).y-td,cs(col,:),'FaceColor',cs(col,:),'FaceAlpha',.3)
                    % get hte ID string for this patch and freq
                    if isfield(school(i).channel(ch),'species')
                        Idstring=[];
                        for sp=1:length(school(i).channel(ch).species)
                            Idstring =[Idstring, ['ID:',school(i).channel(ch).species(sp).speciesID,' fraction:',school(i).channel(ch).species(sp).fraction,';']];
                        end
                    else
                        Idstring='No species ID';% No ID but has frequency so it is ok
                    end
                    if plID
                        text(school(i).x(1),school(i).y(1)-td,Idstring)
                    end
                end
            end
        end
    end
end
% Gavin: can you fix these? I think we need to use a similar structure as
% school and layer, if at all poassible.

% Plot erased regions
% if ~isempty(erased)
%     k = find(ch == [erased.channel.channelID]); % erased data for channel f.
%     if ~isempty(k == 1)
%         for i=1:length(erased.channel(k).x) % loop over each ping with erased samples
%             ping = erased.channel(k).x(i);
%             ranges = erased.channel(k).y{i};
%             for j=1:size(ranges,1) % loop over each contingous block of erased samples
%                 startR = ranges(j,1);
%                 endR = startR + ranges(j,2);
%                 patch([ping ping+1 ping+1 ping], ...
%                     [startR startR endR endR]-td, 'k', ...
%                     'FaceAlpha', 0.8, 'EdgeColor', 'None')
%             end
%         end
%     end
% end
warning('Exclude and erased should work on ping, not time... fix in the reader?')
% if ~isempty(exclude)
%     % Plot exclude regions
%     ym=ylim;
%     maxrange = ym(2);
%     for i=1:length(exclude)
%         [~, startPing] = min(abs(exclude(i).startTime - Sv.pings(f).time));
%         endPing = startPing + exclude(i).numOfPings;
%         patch([startPing startPing endPing endPing], ...
%             [0 maxRange maxRange 0]-td, 'k', 'FaceAlpha', 0.7)
%     end
% end
