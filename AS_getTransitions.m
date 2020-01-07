% calculate state transition frequency for a single trial
% zeke barger 121119 
function result = AS_getTransitions(data, chunk_length, binsize)

binsize = round(binsize/chunk_length); % number of time chunks per bin

result = zeros(9,length(data)/binsize); % results container

% possible pairs
combos = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3];

transitions = zeros(1,(length(data)-1));

for j = 1:9
    transitions(strfind(data, combos(j,:))) = j;
end

startpts = 0:binsize:(length(data)-binsize); startpts(1) = 1;
endpts = (binsize:binsize:length(data)) - 1;

% calculate probabilities in each bin
for i = 1:size(result,2)
    for j = 1:9
        result(j,i) = sum(transitions(startpts(i):endpts(i)) == j) / (endpts(i) - startpts(i) + 1);
    end
end
