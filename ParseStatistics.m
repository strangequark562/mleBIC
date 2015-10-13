function summary = ParseSummary(results,models)

numModels = length(models);
numTracks = length(results);

% intialize variables
popFrac = zeros(1,numModels);
modelIndex = zeros(numTracks,1);
modelProb = zeros(numTracks,numModels);
D = []; sigma = []; v = []; L = []; alpha = []; 
stateD = cell(numModels,1); 
stateSigma = cell(numModels,1);
for i = 1:numModels
    stateD{i} = [];
    stateSigma{i} = [];
end

% get properties
for n = 1:numTracks
    index = find(strcmp(results(n).model,models)==1);
    
    % optimal model index
    modelIndex(n) = index;
    
    % model prob
    for i = 1:numModels
        modelProb(n,i) = results(n).estimates(i).modelProb;
    end
    
    % population fractions
    popFrac(index) = popFrac(index) + 1;

    % optimal estimates
    mu = results(n).mu;
    switch results(n).model
        case 'Immobile'
            stateSigma{1} = [stateSigma{1}; mu];
        case 'Normal'
            stateD{2} = [stateD{2}; mu(1)];
            stateSigma{2} = [stateSigma{2}; mu(2)];
            D = [D; mu(1)];
            sigma = [sigma; mu(2)];
        case 'Driven'
            stateD{3} = [stateD{3}; mu(1)];
            v = [v; sqrt(mu(2)^2+mu(3)^2)];
            stateSigma{3} = [stateSigma{3}; mu(4)];
            D = [D; mu(1)];
            sigma = [sigma; mu(4)];
        case 'Driven'
            stateD{4} = [stateD{4}; mu(1)];
            L = [L; mu(2)];
            stateSigma{4} = [stateSigma{4}; mu(3)];
            D = [D; mu(1)];
            sigma = [sigma; mu(3)];
        case 'Driven'
            stateD{5} = [stateD{5}; mu(1)];
            alpha = [alpha; mu(2)];
            stateSigma{5} = [stateSigma{5}; mu(3)];
            D = [D; mu(1)];
            sigma = [sigma; mu(3)];
        otherwise
    end
end
popFrac = popFrac/length(results);

summary.models = models;
summary.numModels = numModels;
summary.numTracks = numTracks;
summary.popFrac = popFrac;
summary.modelIndex = modelIndex;
summary.modelProb = modelProb;

ensemble.meanD = mean(D);
ensemble.stdD = std(D);
ensemble.meanSigma = mean(sigma);
ensemble.stdSigma = std(sigma);
ensemble.D = D;
ensemble.sigma = sigma;

for i = 1:numModels
    eval(['state.' models{i} '.popFrac = popFrac(i);']);

    if i > 1
        eval(['state.' models{i} '.meanD = mean(stateD{i});']);
        eval(['state.' models{i} '.stdD = std(stateD{i});']);
        eval(['state.' models{i} '.D = stateD{i};']);
    end    
    eval(['state.' models{i} '.meanSigma = mean(stateSigma{i});']);
    eval(['state.' models{i} '.stdSigma = std(stateSigma{i});']);
    eval(['state.' models{i} '.sigma = stateSigma{i};']);

    if i == 3
        eval(['state.' models{i} '.meanV = mean(v);']);
        eval(['state.' models{i} '.stdV = std(v);']);
        eval(['state.' models{i} '.v = v;']);
    elseif i == 4
        eval(['state.' models{i} '.meanL = mean(L);']);
        eval(['state.' models{i} '.stdL = std(L);']);
        eval(['state.' models{i} '.L = L;']);
    elseif i == 5
        eval(['state.' models{i} '.meanAlpha = mean(alpha);']);
        eval(['state.' models{i} '.stdAlpha = std(alpha);']);
        eval(['state.' models{i} '.alpha = alpha;']);
    end
end





