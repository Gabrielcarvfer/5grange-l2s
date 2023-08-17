function [BER, BLER] = runAWGNSimulation(seed, code_rate, block_size, constelation_size, num_blocks, SNR)
% runAWGNSimulation Run the AWGN simulation 
% 
%   seed                - Simulation seed
%   code_rate           - Code rate
%   block_size          - Size of the encoded block
%   constelation_size   - Constelation size
%   num_blocks          - Number of simulated blocks
%   SNR                 - SNR value
%

rng(seed);

%% Load the Code folder
if ~isdeployed
    addpath('./5g-nr-ldpc');
    addpath('./5g-nr-ldpc/codes');
end

% Get the LPDC structure
ldpc_data = ldpcGet(block_size, code_rate);

% Generate the bits
fprintf('Generate data');
tic;
data = randi([0 1], num_blocks, ldpc_data.numInfBits);
fprintf(' - %.2f ms\n', 1000*toc);

% Encoded data
fprintf('Encode data');
rx_data = ldpcEncode(data, ldpc_data);
fprintf(' - %.2f ms\n', 1000*toc);

% Add zero-padding if needed
padding_size = ceil(size(rx_data, 2)/log2(constelation_size))*log2(constelation_size) - size(rx_data, 2);
if padding_size > 0
    rx_data = [rx_data, zeros(size(rx_data, 1), padding_size)];
end

% Modulated signal (baseband)
fprintf('Modulate data');
rx_data = qammod(rx_data.', constelation_size, 'InputType', 'bit', 'UnitAveragePower', true).';
fprintf(' - %.2f ms\n', 1000*toc);

% Add AWGN noise with a given SNR
fprintf('AWGN channel');
rx_data = awgn(rx_data, SNR);
fprintf(' - %.2f ms\n', 1000*toc);

% Approx LLR demapping
fprintf('Demodulate received data');
rx_data = qamdemod(rx_data.', constelation_size, 'OutputType', 'approxllr', 'UnitAveragePower', true).';
fprintf(' - %.2f ms\n', 1000*toc);

% Decoded data
fprintf('Decode received data');
decoded_data = ldpc_decode(rx_data, ldpc_data);
fprintf(' - %.2f ms\n', 1000*toc);

% Calculate BLER and BER
diffs = (decoded_data ~= data);
BLER = sum(sum(diffs, 2) > 0)/num_blocks;
BER = sum(diffs(:)) / (num_blocks * ldpc_data.numInfBits);

% Save information
filename = sprintf('awgn_s_%d_cr_%.2f_bs_%d_qam_%d_snr_%.1f.mat', seed, code_rate, block_size, constelation_size, SNR);
save(fullfile('results', filename), 'seed', 'code_rate', 'block_size', 'constelation_size', ...
    'SNR', 'num_blocks', 'BLER', 'BER', ...
    '-v7.3')

end
