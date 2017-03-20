function status=LSSSreader_comparesnapfiles(file1,file2)
%
% this function takes the output from LSSSreader_readsnapfiles and compares
% the content of the two files.
%
% inout :
% file1 : File to be compared
% file2 : File to be compared
%
% Output
% status.equal : False if the content is different
%

status =1;
% Loop over structures
if length(file1)~=length(file2)
    msg = 'Different number of polygons';
    status = 2;
else
%     for i=1:length(file1)
%         file1{i}
%     end
end