function status=LSSSreader_comparesnapfiles(file1,file2)
%
% this function takes the output structure from LSSSreader_readsnapfiles
% and compares the content of the two files.
%
% inout :
% file1 : Output from LSSSreader_readsnapfiles
% file2 : Output from LSSSreader_readsnapfiles
%
% Output
% status : Status of the comparison
% status = 1 : No changes detected
% status = 2 : Different numbers of regions
% status = 3 : Same number of regions but different area

status = 1;
% Loop over structures
if length(file1)~=length(file2)
    msg = 'Different number of polygons';
    status = 2;
else
    for i=1:length(file1)
        file1{i}
    end
end

