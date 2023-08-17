%% Clean up
clear all;
clc;

%% Init
if ~isdeployed
    addpath('./codes');
end

% Constellation size
M = 256;

% LDPC config
blkSize = 8192;
codeRate = 1/3;

% Get LDPC struct
LDPC = ldpcGet(blkSize, codeRate);

% Simulation parameters
ebno = 80;
numIter = 10;
numErr = 0;

% Convert E_b/N_0 to some SNR
snr = ebno + 10*log10(log2(M)) + 10*log10(codeRate);

%% Simulate
for i = 1:numIter
    
    % Generate random data
    data = randi([0 1], 1, LDPC.numInfBits);

    % Encode
    dataEnc = ldpcEncode(data, LDPC);

    % QAM mapping
    dataMod = qammod(dataEnc(:), M, 'InputType', 'bit', 'UnitAveragePower', true);

    % AWGN
    dataRx = awgn(dataMod, snr);

    % LLR demapping
    dataLlr = qamdemod(dataRx, M, 'OutputType', 'llr', 'UnitAveragePower', true, 'NoiseVariance', 1/10^(snr/10));

    % Decode
    dataHat = ldpcDecode(dataLlr', LDPC);

    % Count number of bit errors
    numErr = numErr + sum(abs(dataHat - data));
    
end

%% BER
ber = numErr / (numIter * LDPC.numInfBits)