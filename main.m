close all
clear all
clc

%% Load the Code folder
if ~isdeployed
    addpath('./5g-nr-ldpc');
    addpath('./5g-nr-ldpc/codes');
end

% % %% Setup scenrio parameters:
% % % central carrier frequency [Hz]
% % fc = 700e6;
% % 
% % % light sped [m/s]
% % c = 3e8;
% % 
% % % wavelength [m]
% % wavelength = c/fc; 
% % 
% % % Total bandwidth [Hz]
% % BW = 23.4e6 / 3;
% % 
% % % Time sampling (is 1 over the total band width). 
% % % Keep in mind that, if you make time_sample < ( 1/BW ) 
% % % the total occupied band width will be larger than BW.
% % time_sample = ( 1/BW );

%% QAM configuration
% Size of the QAM constelation 
constelation_size = 64;

%% LDPC configuration
% Code rate
code_rate = 3/4;

% Block size
block_size = 256;

% Get the LPDC structure
ldpc_data = ldpcGet(block_size, code_rate);

% Bit energy per noise ratio
EbN0 = 20; 

% SNR
SNR = 30%EbN0 + 10*log10(log2(constelation_size)) + 10*log10(code_rate);

% Number of blocks
num_blocks = 1e3;

% Generate the bits
fprintf('Generate data');
tic;
data = randi([0 1], num_blocks, ldpc_data.numInfBits);
fprintf(' - %.2f ms\n', 1000*toc);

% Encoded data
fprintf('Encode data');
encoded_data = ldpcEncode(data, ldpc_data);
fprintf(' - %.2f ms\n', 1000*toc);

padding_size = ceil(size(encoded_data, 2)/log2(constelation_size))*log2(constelation_size) - size(encoded_data, 2);
if padding_size > 0
    encoded_data = [encoded_data, zeros(size(encoded_data, 1), padding_size)];
end

% Modulated signal (baseband)
fprintf('Modulate data');
signal = qammod(encoded_data.', constelation_size, 'InputType', 'bit', 'UnitAveragePower', true).';
fprintf(' - %.2f ms\n', 1000*toc);

% Add AWGN noise with a given SNR
fprintf('AWGN channel');
rec_signal = awgn(signal, SNR, 'measured');
fprintf(' - %.2f ms\n', 1000*toc);

% LLR demapping
fprintf('Demodulate received data');
rec_data = qamdemod(rec_signal.', constelation_size, 'OutputType', 'approxllr', 'UnitAveragePower', true).';
fprintf(' - %.2f ms\n', 1000*toc);

% Decoded data
fprintf('Decode received data');
decoded_data = ldpc_decode(rec_data, ldpc_data);
% decoded_data2 = ldpc_decode(rec_data(:, 1:(size(rec_data, 2) - padding_size)), ldpc_data);
fprintf(' - %.2f ms\n', 1000*toc);
% fprintf('Diff == %d\n', sum(decoded_data(:) ~= decoded_data2(:)));

% decoded_data2 = zeros(num_blocks, ldpc_data.numInfBits);
% for i = 1:num_blocks
%     fprintf('Msg: %d\n', i);
%     decoded_data2(i,:) = ldpcDecode(rec_data(i, :), ldpc_data);
% end
% fprintf('Diff: %d\n', sum(sum(decoded_data(:) ~= decoded_data2(:))))

diffs = (decoded_data ~= data);
BLER = sum(sum(diffs, 2) > 0)/num_blocks;

BER = sum(diffs(:)) / (num_blocks * ldpc_data.numInfBits);

filename = sprintf('awgn_cr_%.2f_bs_%d_qam_%d_snr_%.1f.mat', code_rate, block_size, constelation_size, SNR);
save(fullfile('results', filename), 'code_rate', 'block_size', 'constelation_size', 'SNR', 'num_blocks', 'BLER', 'BER', ...
    'data', 'decoded_data', ...
    '-v7.3')
