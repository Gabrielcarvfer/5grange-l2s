close all
clear all
clc

folder = './results';
seed = 0;

block_size_values = [256, 512, 1024, 2048, 4096];
code_rate_values = [1/3, 1/2, 2/3, 3/4, 5/6];
constelation_size_values = [2, 4, 16, 64, 256];
SNR_values = -10:0.25:30;

BLER = zeros(length(block_size_values), length(code_rate_values), length(constelation_size_values), length(SNR_values));
BER = zeros(length(block_size_values), length(code_rate_values), length(constelation_size_values), length(SNR_values));

block_idx = 1;
for block_size = block_size_values
    fprintf('Processing Block Size = %d\n', block_size);
    code_rate_idx = 1;
    for code_rate = code_rate_values
        fprintf('  Processing Code Rate = %.2f\n', code_rate);
        constelation_idx = 1;
        for constelation_size = constelation_size_values
            fprintf('    Processing modulation = %d-QAM\n', constelation_size);
            snr_idx = 1;
            for SNR = SNR_values
%                fprintf('      Processing SNR = %.2f\n', SNR);
                filename = sprintf('awgn_s_%d_cr_%.2f_bs_%d_qam_%d_snr_%.1f.mat', seed, code_rate, block_size, constelation_size, SNR);
                filepath = fullfile(folder, filename);
                
                if exist(filepath, 'file')
                    file = matfile(filepath, 'Writable', false);
                
                    BLER(block_idx, code_rate_idx, constelation_idx, snr_idx) = file.BLER;
                    BER(block_idx, code_rate_idx, constelation_idx, snr_idx) = file.BER;
                else
                    BLER(block_idx, code_rate_idx, constelation_idx, snr_idx) = NaN;
                    BER(block_idx, code_rate_idx, constelation_idx, snr_idx) = NaN;
                end                    
                
                snr_idx = snr_idx + 1;
            end
            constelation_idx = constelation_idx + 1;
        end
        code_rate_idx = code_rate_idx + 1;
    end
    block_idx = block_idx + 1;
end

save('test.mat', 'BLER', 'BER', 'block_size_values', 'code_rate_values', 'constelation_size_values', 'SNR_values', '-v7.3');

figure;
semilogy(SNR_values, BLER, 'LineWidth', 2)
grid on;
xlabel('SNR (dB)');
ylabel('BLER');

figure;
semilogy(SNR_values, BER, 'LineWidth', 2)
grid on;
xlabel('SNR (dB)');
ylabel('BER');
