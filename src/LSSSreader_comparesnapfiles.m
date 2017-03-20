function status=LSSSreader_comparesnapfiles(region1,region2)
%
% this function takes the output structure from LSSSreader_readsnapfiles
% and compares the content of the two files.
%
% inout :
% region1 : Output from LSSSreader_readsnapfiles
% region2 : Output from LSSSreader_readsnapfiles
%
% Output
% status : Status of the comparison
% status = 1 : No changes detected
% status = 2 : Different number of polygons
% status = 3 : Different areas in at least one polygon (relative difference > 1e-3)
% status = 4 : Different number of polygons AND Different areas in at least one polygon
%

tol = 1e-3;

status1 = 0;
status2 = 0;
msg1 = '';
msg2 = '';

% Loop over structures
if length(region1)~=length(region2)
    msg1 = 'Different number of polygons';
    status1 = 1;
else
    for i=1:length(region1)
        
        A1=polyarea(region1(i).x,region1(i).y);
        A2=polyarea(region2(i).x,region2(i).y);
        if (A1-A2)/A1>tol
            msg2 = 'Different areas in polygon';
            status2 = 1;
        end
    end
end

status = 1 + status1 + status2*2;
