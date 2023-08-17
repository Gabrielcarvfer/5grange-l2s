function channel_per_block = generate_channel_samples(qamOrder, channelModel, fecType, blockSize, codeRate, numerology_id, centerFrequency, mobileSpeed, num_samples, initial_seed)


if ~exist('initial_seed', 'var')
    rng('shuffle')
    initial_seed = randi([0 1e5]);
end

%% Channel coding parameters
if strcmp(fecType, 'LDPC')          % LDPC code settings
    % Get LDPC struct
    LDPC = ldpcGet(blockSize, codeRate);
    actualK = LDPC.numInfBits;
    actualN = LDPC.numTotBits;
    LDPC.iterations = 20;
    LDPC.decType = 'SPA';
elseif strcmp(fecType, 'Polar')     % Polar code settings
    polarN = blockSize;
    polarK = round(polarN * str2num(codeRate));
    polarCodeRate = polarK/polarN;
    polarListSize = 8;
    actualK = polarK;
    actualN = polarN;
else
    error('Error: Unknown fecType = "%s"!', fecType);
end

%% Waveform parameters
% Get all numerologies of FlexNOW
flexNOW = flexnow_getWaveformConfig();

gfdm_params = flexNOW(numerology_id + 1);

% Number of RBs
numRBs = gfdm_params.NRB;

% Get the number of data and pilot REs per RB
[numDataREsPerRB, ~, numPilotsPerRB] = flexnow_getNumRE(flexNOW, numerology_id, 0, 1);

% Obtain the data and pilot indexes for all RBs
data_indexes  = zeros(numDataREsPerRB, numRBs);
pilot_indexes = zeros(numPilotsPerRB,  numRBs);
for rb = 1:numRBs
    [~, data_indexes(:, rb), pilot_indexes(:, rb)] = flexnow_getResurceMapIndex(flexNOW, numerology_id, rb - 1, 1);
end
data_indexes  = data_indexes + 1;
pilot_indexes = pilot_indexes + 1;

numBitsPerSymbol = log2(qamOrder);

% Calculate the size of the zero padding
numSymbolsPerBlock = ceil(actualN / numBitsPerSymbol);

% Get the number of RBs used to transmit a single block
numRBsPerBlock = ceil(numSymbolsPerBlock / numDataREsPerRB);

% Calculate the number of filler bits that shall be added to fill an integer number of RBs
numFillerBits = numRBsPerBlock * numDataREsPerRB * numBitsPerSymbol - actualN;

% Number of blocks
numBlocks = floor(numRBs/numRBsPerBlock);

% Number of RBs fulfilled
numFilledRBs = numBlocks * numRBsPerBlock;

numSubcarriers            = gfdm_params.K;        % Number of Subcarriers
numSubsymbols             = gfdm_params.M;        % Number of subsymbos
cpDuration                = gfdm_params.Ncp;      % Cyclic prefix size
csDuration                = gfdm_params.Ncs;      % Cyclic sufix size
symbolDuration            = numSubcarriers * numSubsymbols;
numGFDMSymbolsPerSubframe = gfdm_params.B;
gfdmSymbolLength          = symbolDuration + cpDuration + csDuration;
pilotFreqSpacing          = gfdm_params.pilotFreqSpacing;
pilotTimeSpacing          = gfdm_params.pilotTimeSpacing;

FFTSize  = symbolDuration;

samplingFrequency = gfdm_params.Fs;
switch (channelModel)
    case {'CDL_A', 'CDL_D'}
        h = FiveGRangeWrapper(channelModel, centerFrequency, mobileSpeed, samplingFrequency, 0, 1);
end

channel_per_block = nan(numBlocks * ceil(num_samples/numBlocks), actualN/numBitsPerSymbol);
seed = initial_seed;

num_channel_blocks = 0;
while num_channel_blocks < num_samples
    switch (channelModel)
        case {'CDL_A', 'CDL_D'}
            % Change the seed of the channel for a new uncorr. channel
            % Deterministically choose seed value
            h.Seed(seed);

            pad = zeros(2 * cpDuration, 1);
            impResp = h.Filter([pad; repmat(circshift([1; zeros(gfdmSymbolLength - 1, 1)], FFTSize/2), numGFDMSymbolsPerSubframe, 1)]);
            impResp = circshift(impResp((numel(pad) + 1):end), -FFTSize/2);
            impResp = reshape(impResp, [], numGFDMSymbolsPerSubframe);
            impResp = impResp(1:FFTSize, :);
        case {'AWGN'}
            impResp = [ones(1, numGFDMSymbolsPerSubframe); zeros(FFTSize - 1, numGFDMSymbolsPerSubframe)];
    end
    actualChannel = fft(impResp);

    channel_sample = actualChannel(data_indexes(:, 1:numFilledRBs)); 
    channel_sample = reshape(channel_sample, (actualN + numFillerBits)/numBitsPerSymbol, numBlocks); 
    channel_sample = channel_sample(1:(actualN/numBitsPerSymbol), :).';
    
    channel_per_block(((seed * numBlocks + 1):((seed + 1) * numBlocks)) - initial_seed * numBlocks, :) = channel_sample;
    num_channel_blocks = num_channel_blocks + numBlocks;
    seed = seed + 1;
end

channel_per_block = channel_per_block(1:num_samples, :);
end