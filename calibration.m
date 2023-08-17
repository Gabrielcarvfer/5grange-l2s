clear all
close all
clc

configure_path();

MCS_table = [ % [Modulation, Code rate]
        {  4, '1/24' };   % MCS 1
        {  4, '1/12' };   % MCS 2
        {  4, '1/8'  };   % MCS 3
        {  4, '1/6'  };   % MCS 4
        {  4, '5/24' };   % MCS 5
        {  4, '7/24' };   % MCS 6
        {  4, '3/8'  };   % MCS 7
        { 16, '5/24' };   % MCS 8
        { 16, '1/4'  };   % MCS 9
        { 16, '7/24' };   % MCS 10
        { 16, '3/8'  };   % MCS 11
        { 16, '11/24'};   % MCS 12
        { 16, '13/24'};   % MCS 13
        { 16, '5/8'  };   % MCS 14
        { 16, '3/4'  };   % MCS 15
        { 16, '5/6'  };   % MCS 16
        { 64, '7/12' };   % MCS 17
        { 64, '2/3'  };   % MCS 18
        { 64, '3/4'  };   % MCS 19
        { 64, '19/24'};   % MCS 20
        { 64, '7/8'  };   % MCS 21
        { 64, '11/12'};   % MCS 22
        {256, '3/4'  };   % MCS 23
        {256, '5/6'  };   % MCS 24
        {256, '7/8'  };   % MCS 25
        {256, '11/12'};   % MCS 26
        {256, '23/24'}    % MCS 27
      ];
  
block_size_values        = [256, 512, 1024, 2048, 4096, 8192];
channel_model_values     = {'CDL_A', 'CDL_D'}; % 'AWGN'
EbN0_values              = -10:0.5:60;
numerology_id_values     = [0, 0, 0, 1,  1,  2,  3,   4,   5];
speed_values             = [0, 3, 7, 0, 15, 30, 60, 120, 240];
perfect_csi              = 1;

center_frequency = 700e6;
num_channel_samples = 100;

lin2dB = @(x) 10*log10(x);
dB2lin = @(x) 10.^(x/10);

awgn_file = matfile(sprintf('results_awgn_est_%d_f700.mat', perfect_csi), 'Writable', false);
chan_file = matfile(sprintf('results_est_%d_f700.mat', perfect_csi), 'Writable', false);

test_beta = 100;

for numerology_idx = 1:length(numerology_id_values)
    numerology_id = numerology_id_values(numerology_idx);
    mobileSpeed   = speed_values(numerology_idx);
    
    for mcs_idx = 25%1:length(MCS_table)
        constelation_size = MCS_table{mcs_idx, 1};
        code_rate         = MCS_table{mcs_idx, 2};
        codeRateStr       = strrep(code_rate,'/','-');

        for bs_idx = 6%1:length(block_size_values)
            block_size = block_size_values(bs_idx);

            for chan_idx = 2%1:length(channel_model_values)
                % Channel model
                channel_model = channel_model_values{chan_idx};
                
                % Channel samples (coded blocks x RE)
                channel = generate_channel_samples(constelation_size, channel_model, 'Polar', block_size, code_rate, numerology_id, center_frequency, mobileSpeed, num_channel_samples);
                
                % SNR and BLER for the AWGN
                awgn_snr_dB = squeeze(awgn_file.SNR(mcs_idx, bs_idx, :));
                awgn_bler   = squeeze(awgn_file.BLER(mcs_idx, bs_idx, :));
                
                % SNR and BLER for the chosen channel model
                chan_snr_dB = squeeze(chan_file.SNR(numerology_idx, mcs_idx, bs_idx, chan_idx, :));
                chan_bler   = squeeze(chan_file.BLER(numerology_idx, mcs_idx, bs_idx, chan_idx, :));
                
                % Calculate BLER of the AWGN for a small SNR step
                snr_fine_dB = -30:0.001:100;
                awgn_bler_fine = interp1(awgn_snr_dB, awgn_bler, snr_fine_dB, 'pchip', 'extrap');
                                
                % Calculate BLER of the channel for a small SNR step
                chan_bler_fine = interp1(chan_snr_dB, chan_bler, snr_fine_dB, 'pchip', 'extrap');
                
                % Calculate all SNRs considering the generated channel
                % models. The set of generated channels are converted into
                % SNR considering one target SNR and stacked into a large
                % matrix
                noise_var  = dB2lin(-awgn_snr_dB);
                snr_values = repmat(abs(channel).^2, length(awgn_snr_dB), 1) ./ repmat(reshape(repmat(noise_var.', size(channel, 1), 1), [], 1), 1, size(channel, 2));
                
                % Calculate the mean SNR per coded block
                mean_snr_per_block_dB = lin2dB(mean(snr_values, 2));
                
                % From the mean SNR, obtain the BLER using curve of the
                % channel model
                chan_bler_calibration = interp1(chan_snr_dB, chan_bler, mean_snr_per_block_dB, 'pchip', 'extrap');
                
                % Select the blocks within the range of acceptable BLER values
                calib_idx = find(chan_bler_calibration >= 0.01 & chan_bler_calibration <= 0.9);
                                
                % Find the SNR values in the AWGN BLER curve with the same
                % BLER value of the channel simulation
                snr_target_calibration_fine_dB = nan(size(chan_bler_calibration, 1), 1);
                for i = 1:length(snr_target_calibration_fine_dB)
                    [~, idx] = min(abs(awgn_bler_fine - chan_bler_calibration(i)).^2);
                    snr_target_calibration_fine_dB(i) = snr_fine_dB(idx);
                end
                
                % Individual SNR values for only the selected blocks
                snr_fine = snr_values(calib_idx,:); 
                snr_target_calibration_fine_dB = snr_target_calibration_fine_dB(calib_idx);
                
                % From the mean SNR, obtain the BLER using the AWGN curve
                awgn_bler_calibration = interp1(awgn_snr_dB, awgn_bler, mean_snr_per_block_dB, 'pchip', 'extrap');
                
                % Find a beta that minimizes the MSE of the effective SNR
                % and the target SNR
                calibrated_beta_snr = fmincon(@(beta) mean(abs(...
                    lin2dB(snr_target_calibration_fine_dB) - arrayfun(@(idx) effective_snr_eesm(snr_fine(idx, :), beta), (1:size(snr_fine, 1)).') ...
                    ).^2), 10, [], [], [], [], 1e-3, 1e5);
%                 calibrated_beta_snr = fmincon(@(beta) mean(abs(...
%                     snr_target_calibration_fine_dB - dB2lin(arrayfun(@(idx) effective_snr_eesm(snr_fine(idx, :), beta), (1:size(snr_fine, 1)).')) ...
%                     ).^2), 10, [], [], [], [], 1e-3, 1e5);
                
                calibrated_beta_bler = fmincon(@(beta) mean(abs( ...
                    log10(chan_bler_calibration(calib_idx)) - max(log10(interp1(snr_fine_dB, awgn_bler_fine, lin2dB(arrayfun(@(idx) effective_snr_eesm(snr_fine(idx, :), beta), (1:size(snr_fine, 1)).')), 'pchip')), -100) ...
                    ).^2), 10, [], [], [], [], 1e-3, 1e5);
                
%                 calibrated_beta_bler = fmincon(@(beta) mean(abs( ...
%                     chan_bler_calibration(calib_idx) - interp1(snr_fine_dB, awgn_bler_fine, lin2dB(arrayfun(@(idx) effective_snr_eesm(snr_fine(idx, :), beta), (1:size(snr_fine, 1)).')), 'pchip') ...
%                     ).^2), 10, [], [], [], [], 1e-3, 1e5);
                
                %[calibrated_ab, errab] = fmincon(@(beta) mean(abs(chan_bler_calibration(calib_idx) - interp1(snr_dB_fine, awgn_bler_fine, lin2dB(arrayfun(@(idx) -beta(1)*log(mean(exp(-snr_fine(idx,:)/beta(2)))), (1:size(snr_fine, 1)))), 'pchip')).^2), [1 1], [], [], [], [], 1e-1, 1e2);
                
                % Mean SNR (blocks x SNR)
                mean_snr_per_block_dB = reshape(mean_snr_per_block_dB, num_channel_samples, length(awgn_snr_dB));
                
                plot_information(snr_fine_dB, awgn_bler_fine, chan_bler_fine, snr_fine, chan_bler_calibration(calib_idx), snr_target_calibration_fine_dB, calibrated_beta_snr, 'SNR calib.');
                plot_information(snr_fine_dB, awgn_bler_fine, chan_bler_fine, snr_fine, chan_bler_calibration(calib_idx), snr_target_calibration_fine_dB, calibrated_beta_bler, 'BLER calib.');
                plot_information(snr_fine_dB, awgn_bler_fine, chan_bler_fine, snr_fine, chan_bler_calibration(calib_idx), snr_target_calibration_fine_dB, test_beta, 'Test calib.');
                
                
                
%                 % Calculate the effective SNR for the calibrated beta and a
%                 % test beta
%                 effective_snr_dB = nan(num_channel_samples, length(awgn_snr_dB));
%                 snr_test_beta_dB = effective_snr_dB;
%                 
%                 for snr_idx = 1:length(awgn_snr_dB)
%                     snr_dB = (abs(channel).^2)./noise_var(snr_idx);
%                     
%                     effective_snr_dB(:, snr_idx) = arrayfun(@(idx) lin2dB(effective_snr_eesm(snr_dB(idx, :), calibrated_beta)), (1:size(snr_dB, 1)).');
%                     snr_test_beta_dB(:, snr_idx) = arrayfun(@(idx) lin2dB(effective_snr_eesm(snr_dB(idx, :), test_beta)), (1:size(snr_dB, 1)).');
%                 end

                
                return;
                for snr_idx = 1:length(awgn_snr_dB) 
                    bler_axis = interp1(awgn_snr_dB, awgn_bler, mean_snr_per_block_dB(:, snr_idx), 'pchip', 'extrap');
                    
                    semilogy(effective_snr_dB(:, snr_idx), bler_axis, '.b','DisplayName', sprintf('Calibrated \\beta = %.3f', calibrated_beta));
%                     semilogy(effective_snr_dB{snr_idx}, awgn_bler(snr_idx) * ones(size(effective_snr_dB{snr_idx})), '.b','DisplayName', sprintf('Calibrated \\beta = %.3f', calibrated_beta));
%                     semilogy(effective_snr_dB{snr_idx}, bler_awgn_effective_snr, '.b','DisplayName', sprintf('Calibrated \\beta = %.3f', calibrated_beta));
                    hold on;
                    semilogy(snr_test_beta_dB(:, snr_idx), bler_axis, '.r', 'DisplayName', sprintf('\\beta = %g', test_beta));
%                     semilogy(snr_test_beta_dB{snr_idx}, awgn_bler(snr_idx) * ones(size(effective_snr_dB{snr_idx})), '.r', 'DisplayName', sprintf('\\beta = %g', test_beta));
%                     semilogy(snr_test_beta_dB{snr_idx}, bler_awgn_test_beta, '.r', 'DisplayName', sprintf('\\beta = %g', test_beta));
                end
                semilogy(awgn_snr_dB, awgn_bler, 'k', 'LineWidth', 2, 'DisplayName', 'AWGN');
                semilogy(chan_snr_dB, chan_bler, '--magenta', 'LineWidth', 2, 'DisplayName', 'Channel');
                grid on;
                xlabel('SNR (dB)');
                ylabel('BLER');
                ylim([1e-4 1]);
                title(sprintf('\\beta calibration for MCS: %d - Block: %d - Channel: %s', mcs_idx, block_size, channel_model));
                curves = get(gca, 'Children');
                legend(curves(1:4), 'Location', 'Best');
                
                return;
            end
        end
    end
end
