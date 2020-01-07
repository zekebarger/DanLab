% zeke barger 121519
% check the files in a dirlist to see if there will be any problems during
% further analysis
% input
% dirlist: the list of data folders to check
% varNames: cell array of variable names to look for
% output
% problemString: describes any problems, if they are found
%   empty if everything looks good
% fileNames: list (in order of dirlist entries) of lists (each in order of 
%   varNames) of paths to files 

function [problemString, fileNames] = AS_checkList(dirlist, varNames)
nEntries = length(dirlist);
% default values
problemString = {};
fileNames = cell(1, nEntries);

% if it's empty, that's fine 
if nEntries == 0
    return
end

% otherwise, look thru each entry
for i = 1:nEntries
    % collect information about the contents of the dirlist
    [entryProblemString, entryFileNames] = AS_checkEntry(dirlist{i}, varNames);
    % append to our problem string
    problemString = [problemString; entryProblemString];
    % add list of filenames to the list
    fileNames{i} = entryFileNames;
end