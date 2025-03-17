% dft_resolution_binwidth.m - Simulating Lower Resolution with and without Zero-Padding
clc; clear; close all;

%% Parameters
sampling_rate = 128;  % Fixed Sampling Rate
FFT_SIZE = 128;       % Fixed FFT Size for Zero-Padding Cases

% Different input sample sizes to compare (Lower N -> Worse Resolution)
sample_sizes = [128, 32, 16, 8];

% Create figure
fig = figure('Position', [100, 100, 1200, 800]);
fig.WindowState = 'maximized';
for idx = 1:length(sample_sizes)
    N = sample_sizes(idx);  % Current number of samples
    t = (0:N-1) / sampling_rate;  % Time vector
    signal = sin(2 * pi * 5 * t) + sin(2 * pi * 6 * t);  % Signal

    % --- CASE 1: No Zero-Padding (FFT Size = N) ---
    X_no_pad = fft(signal, N);
    X_no_pad = X_no_pad / sqrt(N);
    magnitude_no_pad = abs(X_no_pad);
    frequencies_no_pad = (0:N-1) * (sampling_rate / N);

    % --- CASE 2: Zero-Padding (FFT Size = Fixed 128) ---
    signal_padded = [signal, zeros(1, FFT_SIZE - N)];
    X_pad = fft(signal_padded, FFT_SIZE);
    X_pad = X_pad / sqrt(FFT_SIZE);
    magnitude_pad = abs(X_pad);
    frequencies_pad = (0:FFT_SIZE-1) * (sampling_rate / FFT_SIZE);

    % --- Plot signal ---
    subplot(4, 2, 2*idx - 1);
    scatter(t*sampling_rate,signal,'filled');
    xlim([0 FFT_SIZE])
    title(sprintf('Signal with N = %d', N));
    grid on;
    set(gca,'FontSize',12) %and other properties
    xlabel('Sample number','FontSize',10);
    ylabel('Amplitude','FontSize',10);

    % --- Plot No Zero-Padding (Lower FFT Resolution) ---
    subplot(4, 2, 2*idx);
    stem(frequencies_no_pad, magnitude_no_pad, 'LineWidth', 1.5);
    title(sprintf('DFT with N = %d', N));
    grid on;
    xline(sampling_rate/2, '--', 'Nyquist Limit','LabelVerticalAlignment', 'middle');
    set(gca,'FontSize',12) %and other properties
    xlim([0 140])
    xlabel('Frequency (Hz)','FontSize',10);
    ylabel('Magnitude','FontSize',10);

end

print(fig, 'dft_resolution_binwidth_signal.png', '-dpng', '-r150'); % 600 DPI high resolution

% Create figure
fig = figure('Position', [100, 100, 1200, 800]);
fig.WindowState = 'maximized';
for idx = 1:length(sample_sizes)
    N = sample_sizes(idx);  % Current number of samples
    t = (0:N-1) / sampling_rate;  % Time vector
    signal = sin(2 * pi * 5 * t) + sin(2 * pi * 6 * t);  % Signal

    % --- CASE 1: No Zero-Padding (FFT Size = N) ---
    X_no_pad = fft(signal, N);
    X_no_pad = X_no_pad / sqrt(N);
    magnitude_no_pad = abs(X_no_pad);
    frequencies_no_pad = (0:N-1) * (sampling_rate / N);

    % --- CASE 2: Zero-Padding (FFT Size = Fixed 128) ---
    signal_padded = [signal, zeros(1, FFT_SIZE - N)];
    X_pad = fft(signal_padded, FFT_SIZE);
    X_pad = X_pad / sqrt(FFT_SIZE);
    magnitude_pad = abs(X_pad);
    frequencies_pad = (0:FFT_SIZE-1) * (sampling_rate / FFT_SIZE);

    % --- Plot No Zero-Padding (Lower FFT Resolution) ---
    subplot(4, 2, 2*idx - 1);
    stem(frequencies_no_pad, magnitude_no_pad, 'LineWidth', 1.5);
    title(sprintf('DFT with N = %d (No Zero-Padding)', N));
    grid on;
    xline(sampling_rate/2, '--', 'Nyquist Limit','LabelVerticalAlignment', 'middle');
    set(gca,'FontSize',12) %and other properties
    xlim([0 140])
    xlabel('Frequency (Hz)','FontSize',10);
    ylabel('Magnitude','FontSize',10);

    % --- Plot with Zero-Padding (Higher FFT Resolution) ---
    subplot(4, 2, 2*idx);
    stem(frequencies_pad, magnitude_pad, 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    title(sprintf('DFT with N = %d (Zero-Padded to %d)', N, FFT_SIZE));
    grid on;
    xline(sampling_rate/2, '--', 'Nyquist Limit', 'LabelVerticalAlignment', 'middle');
    set(gca,'FontSize',12) %and other properties
    xlim([0 140])
    xlabel('Frequency (Hz)','FontSize',10);
    ylabel('Magnitude','FontSize',10);

end


% Save as high-resolution PNG
%saveas(fig, 'dft_resolution_binwidth.png','-dpng','-r600');
print(fig, 'dft_resolution_binwidth.png', '-dpng', '-r150'); % 600 DPI high resolution
%disp('High-resolution figure saved as dft_resolution_binwidth.png');

%fig2plotly(gcf)

