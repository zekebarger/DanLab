% zeke barger 121519
% look through a folder to see if the .mat files inside collectively
%   contain exactly one variable with each name
% input
% entry: name of the folder to search
% varNames: cell array of variable names to look for
% output
% problemString: empty if there's no problems, otherwise describes the
%   problem
% fileNames: cell array of paths to files containing the variables, in the
%   same order as varNames

function [problemString, fileNames] = AS_checkEntry(entry, varNames)
nVars = length(varNames);

% default values
problemString = {};
fileNames = cell(1, nVars);

% if it's empty, that's fine 
if isempty(entry) || nVars == 0
    return
end

% otherwise, look thru folders in the entry 
% if there's not exactly one varName variable in the
% folder, that's a problem

% get contents of folder
d = struct2cell(dir(entry));
% see which files are .mat files
mat = contains(d(1,:),'.mat');
% count .mat files in the folder
nFiles = sum(mat);
% if there aren't any files
if nFiles == 0
    problemString{1,1} = [entry, ' is empty'];
    return
end

% get indices of .mat files
matIdx = find(mat);
% keep track of which variables have already been found
varsFound = zeros(1,nVars);
% for each file
for i = 1:nFiles
    % get full name of the file
    fName = [d{2,matIdx(i)},'\',d{1,matIdx(i)}];
    % for each variable
    for j = 1:nVars
        % if variable is in the file
        if ~isempty(whos('-file',fName, varNames{j}))
            % note that this variable has been found
            varsFound(j) = varsFound(j)+1;
            % store location of this variable
            fileNames{j} = fName;
            % if we already found it though...
            if varsFound(j) == 2
                problemString{end+1,1} = [entry, ' contains more than one ',varNames{j},' variable'];
            end
        end
    end
end

% if we didn't find all the variables
if ~all(varsFound)
    % get the variables we didn't find
    notFound = find(~varsFound);
    for i = 1:length(notFound)
        problemString{end+1,1} = [entry, ' is missing a ',varNames{notFound(i)},' variable'];
    end
end