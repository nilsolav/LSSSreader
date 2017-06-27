function [school,layer,exclude,erased]=LSSSreader_readsnapfiles(file)
% Reads the LSSS nap and work files and generates polygons for each region
% and school
%
% layerInterpretation=LSSSreader_readsnapfiles(file)
%
% Input:
% file : The data file


%         <layer id="0" hasBeenVisisted="true">
%             <restSpecies ID="6"/>
%             <speciesInterpretationRoot>
%                <speciesInterpretationRep frequency="38">
%                   <species ID="6" fraction="1.0"/>
%                </speciesInterpretationRep>
%                <speciesInterpretationRep frequency="70">
%                   <species ID="6" fraction="1.0"/>
%                </speciesInterpretationRep>
%                <speciesInterpretationRep frequency="120">
%                   <species ID="6" fraction="1.0"/>
%                </speciesInterpretationRep>
%                <speciesInterpretationRep frequency="200">
%                   <species ID="6" fraction="1.0"/>
%                </speciesInterpretationRep>
%                <speciesInterpretationRep frequency="364">
%                   <species ID="6" fraction="1.0"/>
%                </speciesInterpretationRep>
%             </speciesInterpretationRoot>
%             <boundaries>

% Old version
%       <schoolRep referenceTime="1.132294810068E9" hasBeenVisisted="true" objectNumber="106">
%          <restSpecies ID="21"/>
%          <speciesInterpretationRoot>
%             <speciesInterpretationRep frequency="38">
%                <species ID="21" fraction="1.0"/>
%             </speciesInterpretationRep>
%          </speciesInterpretationRoot>
%          <boundaryPoints>1834 188.6735 1834 187.71759 1833 186.87169 1832 181.94098 1831 181.55861
% 1830 181.26744 1829 181.25745 1828 181.09627 1827 181.05626 1826 181.89217
% 1825 181.87334 1825 183.21161 1826 183.6128 1827 183.92397 1827 184.49751
% 1826 184.5687 1825 184.35869 1824 184.67107 1824 186.00932 1825 186.46169
% 1826 186.86287 1827 186.98286 1828 187.02287 1829 187.18405 1830 187.95877
% 1831 189.01466 1832 189.20584 1832 185.95576 1833 189.16585 </boundaryPoints>
%       </schoolRep>
%
%
% And we need to suport the new one:
%
%    <schoolInterpretation>
%       <schoolMaskRep referenceTime="1.462268115355E9" hasBeenVisisted="true" objectNumber="6079">
%          <speciesInterpretationRoot>
%             <speciesInterpretationRep frequency="18">
%                <species ID="6010" fraction="1"/>
%             </speciesInterpretationRep>
%             <speciesInterpretationRep frequency="38">
%                <species ID="6010" fraction="1"/>
%             </speciesInterpretationRep>
%             <speciesInterpretationRep frequency="70">
%                <species ID="6010" fraction="1"/>
%             </speciesInterpretationRep>
%             <speciesInterpretationRep frequency="120">
%                <species ID="6010" fraction="1"/>
%             </speciesInterpretationRep>
%             <speciesInterpretationRep frequency="200">
%                <species ID="6010" fraction="1"/>
%             </speciesInterpretationRep>
%             <speciesInterpretationRep frequency="333">
%                <species ID="6010" fraction="1"/>
%             </speciesInterpretationRep>
%          </speciesInterpretationRoot>
%          <pingMask relativePingNumber="205">45.21126 50.080643</pingMask>
%          <pingMask relativePingNumber="206">45.21126 50.080643</pingMask>
%          <pingMask relativePingNumber="207">45.21126 50.080643</pingMask>
%          <pingMask relativePingNumber="208">45.21126 50.080643</pingMask>
%          <pingMask relativePingNumber="209">45.21126 50.080643</pingMask>
%          <pingMask relativePingNumber="210">45.21126 50.080643</pingMask>
%          <pingMask relativePingNumber="211">45.21126 50.080643</pingMask>
%          <pingMask relativePingNumber="212">45.21126 50.080643</pingMask>
%       </schoolMaskRep>

%% Import the snap file
D.snap = xml2struct(file);

% Fix a problem where xml2struct doesn't return a cell array if only one
% item is present. Sort of a hack...
if length(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer) == 1
    l = D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer;
    D.snap.regionInterpretation.layerInterpretation.layerDefinitions = rmfield(D.snap.regionInterpretation.layerInterpretation.layerDefinitions, 'layer');
    D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{1} = l;
end

if isfield(D.snap.regionInterpretation, 'masking')
    if isfield(D.snap.regionInterpretation.masking, 'mask')
        if length(D.snap.regionInterpretation.masking.mask) == 1
            m = D.snap.regionInterpretation.masking.mask;
            D.snap.regionInterpretation.masking = rmfield(D.snap.regionInterpretation.masking, 'mask');
            D.snap.regionInterpretation.masking.mask{1} = m;
        end
    end
end

if isfield(D.snap.regionInterpretation.exclusionRanges, 'timeRange')
    if length(D.snap.regionInterpretation.exclusionRanges.timeRange) == 1
        t = D.snap.regionInterpretation.exclusionRanges.timeRange;
        D.snap.regionInterpretation.exclusionRanges = rmfield(D.snap.regionInterpretation.exclusionRanges, 'timeRange');
        D.snap.regionInterpretation.exclusionRanges.timeRange{1} = t;
    end
end
%% Get the excluded ping ranges
exclude = struct([]);
if isfield(D.snap.regionInterpretation.exclusionRanges, 'timeRange')
    timeRange = D.snap.regionInterpretation.exclusionRanges.timeRange;
    nsE = length(timeRange);
    for i = 1:nsE
        exclude(i).startTime = unixTimeToMatlab(str2double(timeRange{i}.Attributes.start)); %#ok<*AGROW>
        exclude(i).numOfPings = str2double(timeRange{i}.Attributes.numberOfPings);
    end
end

%% Get the erased masks
% Generates a structure called 'erased'. For each named ping, there is one
% or more depth ranges that are to be erased.
erased = struct([]);
if isfield(D.snap.regionInterpretation, 'masking')
    if isfield(D.snap.regionInterpretation.masking, 'mask')
        erased(1).referenceTime = D.snap.regionInterpretation.masking.Attributes.referenceTime;
        nsM = length(D.snap.regionInterpretation.masking.mask); % one for each channel
        for i = 1:nsM
            m = D.snap.regionInterpretation.masking.mask{i};
            erased(1).channel(i).channelID = str2double(m.Attributes.channelID);
            for j = 1:length(m.ping)
                erased(1).channel(i).x(j) = str2double(m.ping{j}.Attributes.pingOffset);
                ranges = str2num(m.ping{j}.Text);
                erased(1).channel(i).y{j} = reshape(ranges, 2, [])';
            end
        end
    end
end
%% Get the schoolInterpretation
school = struct([]);
if isfield(D.snap.regionInterpretation.schoolInterpretation, 'schoolRep')
    nsI= length(D.snap.regionInterpretation.schoolInterpretation.schoolRep);
    for i=1:nsI
        % The boundary
        %        D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep.species.Attributes.ID
        if isfield(D.snap.regionInterpretation.schoolInterpretation.schoolRep{i},'speciesInterpretationRoot') && isfield(D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot,'speciesInterpretationRep')
            % Need to implement more than one region for this one...
            nfreq=length(D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep);
            if nfreq==1
                school(i).channel(1).speciesID =  D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep.species.Attributes.ID;
                school(i).channel(1).fraction =  D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep.species.Attributes.fraction;
                school(i).channel(1).frequency = D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep.Attributes.frequency;
            else
               for fr = 1:nfreq 
                   school(i).channel(fr).speciesID =  D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep{fr}.species.Attributes.ID;
                   school(i).channel(fr).fraction =  D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep{fr}.species.Attributes.fraction;
                   school(i).channel(fr).frequency = D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.speciesInterpretationRoot.speciesInterpretationRep{fr}.Attributes.frequency;
               end
            end
        end
        T=D.snap.regionInterpretation.schoolInterpretation.schoolRep{i}.boundaryPoints.Text;
        dum = str2num(strrep(T,newline,' '));
        school(i).x = dum(1:2:end-1);
        school(i).y = dum(2:2:end);
        school(i).regiontype = 'school';
    end
end

%% Get the layerInterpretation.boundaries.verticalBoundary(ies)
nvB= length(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary);

%verticalBoundary = NaN(nvB,6);
for i = 1:nvB
    % The attributes
    layerInterpretation.boundaries.verticalBoundary(i).referenceTime = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary{i}.verticalBoundaryRep.Attributes.referenceTime);
    layerInterpretation.boundaries.verticalBoundary(i).pingOffset = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary{i}.verticalBoundaryRep.Attributes.pingOffset);
    layerInterpretation.boundaries.verticalBoundary(i).startDepth = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary{i}.verticalBoundaryRep.Attributes.startDepth);
    layerInterpretation.boundaries.verticalBoundary(i).endDepth   = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary{i}.verticalBoundaryRep.Attributes.endDepth);
    % The connections
    layerInterpretation.boundaries.verticalBoundary(i).endConnector   = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary{i}.Attributes.endConnector);
    layerInterpretation.boundaries.verticalBoundary(i).id             = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary{i}.Attributes.id);
    layerInterpretation.boundaries.verticalBoundary(i).startConnector = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.verticalBoundary{i}.Attributes.startConnector);
end

% Get the layerInterpretation.boundaries.curveBoundary(ies)
ncB= length(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary);
for i = 1:ncB
    % The attributes
    layerInterpretation.boundaries.curveBoundary(i).endConnector   = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary{i}.Attributes.endConnector);
    layerInterpretation.boundaries.curveBoundary(i).id             = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary{i}.Attributes.id);
    layerInterpretation.boundaries.curveBoundary(i).startConnector = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary{i}.Attributes.startConnector);
    layerInterpretation.boundaries.curveBoundary(i).referenceTime  = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary{i}.curveRep.pingRange.Attributes.referenceTime);
    layerInterpretation.boundaries.curveBoundary(i).startOffset    = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary{i}.curveRep.pingRange.Attributes.startOffset);
    layerInterpretation.boundaries.curveBoundary(i).numberOfPings  = str2double(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary{i}.curveRep.pingRange.Attributes.numberOfPings);
    layerInterpretation.boundaries.curveBoundary(i).depths = str2num(strrep(D.snap.regionInterpretation.layerInterpretation.boundaries.curveBoundary{i}.curveRep.depths.Text,newline,' ')); %#ok<*ST2NM>
end

% Get the layerInterpretation.connector (s)
nc= length(D.snap.regionInterpretation.layerInterpretation.connectors.connector);
for i = 1:nc
    layerInterpretation.connector(i).depth   = str2double(D.snap.regionInterpretation.layerInterpretation.connectors.connector{i}.connectorRep.Attributes.depth);
    layerInterpretation.connector(i).pingOffset   = str2double(D.snap.regionInterpretation.layerInterpretation.connectors.connector{i}.connectorRep.Attributes.pingOffset);
    layerInterpretation.connector(i).referenceTime   = str2double(D.snap.regionInterpretation.layerInterpretation.connectors.connector{i}.connectorRep.Attributes.referenceTime);
    layerInterpretation.connector(i).id   = str2double(D.snap.regionInterpretation.layerInterpretation.connectors.connector{i}.Attributes.id);
end


% layerInterpretation.layerDefinitions.layer
nL = length(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer);
for i = 1:nL
    layerInterpretation.layer(i).id   = str2double(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.Attributes.id);
    layerInterpretation.layer(i).hasBeenVisisted = D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.Attributes.hasBeenVisisted;
    layerInterpretation.layer(i).restSpecies = 0.0; % IS THIS THE RIGHT DEFAULT VALUE when restSpecies is not in the XML????
    if isfield(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}, 'restSpecies')
        layerInterpretation.layer(i).restSpecies = str2double(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.restSpecies);
    end
    
    
    for j=1:length(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.boundaries.curveBoundary)
        layerInterpretation.layer(i).curveBoundary(j) = str2double(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.boundaries.curveBoundary{j}.Attributes.id);
    end
    for k=1:length(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.boundaries.verticalBoundary)
        layerInterpretation.layer(i).verticalBoundary(k) = str2double(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.boundaries.verticalBoundary{k}.Attributes.id);
    end
    for k=1:length(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.connectors.id)
        layerInterpretation.layer(i).connectors(k) = str2double(D.snap.regionInterpretation.layerInterpretation.layerDefinitions.layer{i}.connectors.id{k}.Attributes.number);
    end
end

%% Debugging plot section
if false
    figure
    hold on
    % Plot connectors
    for j=1:length(layerInterpretation.connector)
        plot(layerInterpretation.connector(j).pingOffset,layerInterpretation.connector(j).depth,'*r')
        text(layerInterpretation.connector(j).pingOffset,layerInterpretation.connector(j).depth,num2str(layerInterpretation.connector(j).id))
    end
    % Plot curved boundaries
    for j=1:length(layerInterpretation.boundaries.curveBoundary)
        t=(1:layerInterpretation.boundaries.curveBoundary(j).numberOfPings)  +layerInterpretation.boundaries.curveBoundary(j).startOffset;
        plot(t,layerInterpretation.boundaries.curveBoundary(j).depths,'k-')
    end
    % Plot vertical boundaries
    for j=1:length(layerInterpretation.boundaries.verticalBoundary)
        t=layerInterpretation.boundaries.verticalBoundary(j).pingOffset;
        cy = [layerInterpretation.boundaries.verticalBoundary(j).startDepth layerInterpretation.boundaries.verticalBoundary(j).endDepth];
        plot([t t],cy)
    end
end

%% Generate regions as matlab polygons

% Loop over layers
for i= 1:length(layerInterpretation.layer)
    % Loop over curved boundaries (and check if they are part of this
    % layer) and get the connections. The connections variable is the key
    % for stitching the boundaries together:
    % The first dimension is the different boundaries for this layer, the
    % second dimension is:
    % connections =[curve=1/vertical=2 j=layerindex_in_matlab_structure
    %               boundary_id start_connector_id end_connector_id]
    
    connections =[];
    for j=1:length(layerInterpretation.boundaries.curveBoundary)
        if ismember(layerInterpretation.boundaries.curveBoundary(j).id,layerInterpretation.layer(i).curveBoundary)
            connections = [connections;[1 j layerInterpretation.boundaries.curveBoundary(j).id layerInterpretation.boundaries.curveBoundary(j).startConnector layerInterpretation.boundaries.curveBoundary(j).endConnector]];
        end
    end
    % Loop over vertical boundaries (and check if they are part of this
    % layer) and get the connections
    for j=1:length(layerInterpretation.boundaries.verticalBoundary)
        if ismember(layerInterpretation.boundaries.verticalBoundary(j).id,layerInterpretation.layer(i).verticalBoundary)
            connections = [connections;[2 j layerInterpretation.boundaries.verticalBoundary(j).id layerInterpretation.boundaries.verticalBoundary(j).startConnector layerInterpretation.boundaries.verticalBoundary(j).endConnector]];
        end
    end
    % It seems that the boundary lines do not form a closed circle. Rather
    % the boundaries forks from a single starting point. I assume that the
    % starting point is the point where there are two equal starting nodes
    % in "connections", i.e. connections(:,3). The following logic peaces
    % this together and forms the polygon for the layer by starting at a
    % node, and then searching for the next, appending it and then, move
    % cycle through the remaing ones.
    % The structure "connections_sorted" is used to order the boundaries
    % such that they form a closed polygon. The first dimension is the
    % boudnaries (similar to the "connections" variable, and the second
    % dimension is :
    % connections_sorted =[curve=1/vertical=2 j=layerindex_in_matlab_structure
    %               boundary_id start_connector_id end_connector_id forward]
    % The logic is that the endnode in the first is the either the
    % startnode (forward==true) or the endnode (forward==false) in the next
    nextnode=1; % Let us start on the first boundary (could be any)
    exind = 1:size(connections,1);
    forward=true;
    connections_sorted=[];
    for j=1:size(connections,1)
        % Find the startingpoint for the boundary
        if forward
            startnode = connections(nextnode,4);
            endnode = connections(nextnode,5);
        else
            startnode = connections(nextnode,5);
            endnode = connections(nextnode,4);
        end
        % Append the connections_sorted variable
        connections_sorted =[connections_sorted; [connections(nextnode,1:3) startnode endnode forward]];
        % Remove the "used" nodes
        exind(nextnode==exind)=[];
        % Find the next node in the connections while removeing the "used
        % nodes"
        residualconnections=connections(exind,4:5);
        ind =  residualconnections == endnode;
        % change "direction"? Needs to change direction
        % multiple times.
        % Go backwards
        if sum(ind(:,2))>0
            forward = false;
        elseif sum(ind(:,1))>0 % Turn forward again
            forward = true;
        end
        nextnode = find(ind(:,1)|ind(:,2));
        nextnode=exind(nextnode);
    end
    
    % Based on the "connections_sorted" variable the x and y pairs are formed
    % for the layer.
    layer(i).x = [];
    layer(i).y = [];
    for j=1:size(connections_sorted,1)
        % The matlab ID for the boundary
        matid = connections_sorted(j,2);
        % If it is a vertical boundary
        if connections_sorted(j,1)==2
            x = [layerInterpretation.boundaries.verticalBoundary(matid).pingOffset layerInterpretation.boundaries.verticalBoundary(matid).pingOffset];
            y = [layerInterpretation.boundaries.verticalBoundary(matid).startDepth layerInterpretation.boundaries.verticalBoundary(matid).endDepth];
        elseif connections_sorted(j,1)==1
            % if it is a curved layer
            x = (1:layerInterpretation.boundaries.curveBoundary(matid).numberOfPings) + layerInterpretation.boundaries.curveBoundary(matid).startOffset;
            y = layerInterpretation.boundaries.curveBoundary(matid).depths;
        else
            error('No layer coordinates')
        end
        
        % Does the layer start with the endnode? Then we need to turn the
        % direction of the points
        if connections_sorted(j,6)~=1
            x=x(end:-1:1);
            y=y(end:-1:1);
        end
        % Append to the previous layer
        layer(i).x = [layer(i).x x];
        layer(i).y = [layer(i).y y];
    end
    % Add information about the layer
    layer(i).fraction  = NaN;
    layer(i).speciesID = NaN;
    layer(i).regiontype = 'region';
end

end

function t = unixTimeToMatlab(tt)
%  convert Java time numer to MATLAB serial time
t = tt / 86400 + datenum(1970, 1, 1);
end

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

if (nargin < 1)
    clc;
    help xml2struct
    return
end

if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
    % input is a java xml object
    xDoc = file;
else
    %check for existance
    if (exist(file,'file') == 0)
        %Perhaps the xml extension was omitted from the file name. Add the
        %extension and try again.
        if (~contains(file,'.xml'))
            file = [file '.xml'];
        end
        
        if (exist(file,'file') == 0)
            error(['The file ' file ' could not be found']);
        end
    end
    %read the xml file
    xDoc = xmlread(file);
end

%parse xDoc into a MATLAB structure
s = parseChildNodes(xDoc);

end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct; textflag = 'Text';
if hasChildNodes(theNode)
    childNodes = getChildNodes(theNode);
    numChildNodes = getLength(childNodes);
    
    for count = 1:numChildNodes
        theChild = item(childNodes,count-1);
        [text,name,attr,childs,textflag] = getNodeData(theChild);
        
        if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
            %XML allows the same elements to be defined multiple times,
            %put each in a different cell
            if (isfield(children,name))
                if (~iscell(children.(name)))
                    %put existsing element into cell format
                    children.(name) = {children.(name)};
                end
                index = length(children.(name))+1;
                %add new element
                children.(name){index} = childs;
                if(~isempty(fieldnames(text)))
                    children.(name){index} = text;
                end
                if(~isempty(attr))
                    children.(name){index}.('Attributes') = attr;
                end
            else
                %add previously unknown (new) element to the structure
                children.(name) = childs;
                if(~isempty(text) && ~isempty(fieldnames(text)))
                    children.(name) = text;
                end
                if(~isempty(attr))
                    children.(name).('Attributes') = attr;
                end
            end
        else
            ptextflag = 'Text';
            if (strcmp(name, '#cdata_dash_section'))
                ptextflag = 'CDATA';
            elseif (strcmp(name, '#comment'))
                ptextflag = 'Comment';
            end
            
            %this is the text in an element (i.e., the parentNode)
            if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                    ptext.(ptextflag) = text.(textflag);
                else
                    %what to do when element data is as follows:
                    %<element>Text <!--Comment--> More text</element>
                    
                    %put the text in different cells:
                    % if (~iscell(ptext)) ptext = {ptext}; end
                    % ptext{length(ptext)+1} = text;
                    
                    %just append the text
                    ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                end
            end
        end
        
    end
end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr)))
    attr = [];
end

%parse child nodes
[childs,text,textflag] = parseChildNodes(theNode);

if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
    %get the data of any childless nodes
    % faster than if any(strcmp(methods(theNode), 'getData'))
    % no need to try-catch (?)
    % faster than text = char(getData(theNode));
    text.(textflag) = toCharArray(getTextContent(theNode))';
end

end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = struct;
if hasAttributes(theNode)
    theAttributes = getAttributes(theNode);
    numAttributes = getLength(theAttributes);
    
    for count = 1:numAttributes
        %attrib = item(theAttributes,count-1);
        %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
        %attributes.(attr_name) = char(getValue(attrib));
        
        %Suggestion of Adrian Wanner
        str = toCharArray(toString(item(theAttributes,count-1)))';
        k = strfind(str,'=');
        attr_name = str(1:(k(1)-1));
        attr_name = strrep(attr_name, '-', '_dash_');
        attr_name = strrep(attr_name, ':', '_colon_');
        attr_name = strrep(attr_name, '.', '_dot_');
        attributes.(attr_name) = str((k(1)+2):(end-1));
    end
end
end