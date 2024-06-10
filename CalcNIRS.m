%% Hadar Granot 208749176 and Maya Wertsman 208915355

function [dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)

    % Validating the inputs
    assert(ischar(dataFile) && isfile(dataFile), 'Data file not found or invalid filename');
    assert(ischar(extinctionCoefficientsFile) && isfile(extinctionCoefficientsFile), 'Extinction coefficients file not found or invalid filename');
    assert(ischar(DPFperTissueFile) && isfile(DPFperTissueFile), 'DPF per tissue file not found or invalid filename');
    assert(ischar(relDPFfile) && isfile(relDPFfile), 'Relative DPF file not found or invalid filename');
    assert(isnumeric(SDS) && isscalar(SDS) && SDS > 0, 'SDS must be a positive scalar');
    assert(ischar(tissueType) && ~isempty(tissueType), 'Tissue type must be a non-empty string');
    assert(isnumeric(plotChannelIdx) && all(plotChannelIdx > 0) && isvector(plotChannelIdx), 'Plot channel indices must be a vector of positive integers');
    
    % Loading the data
    data = load(dataFile);
    requiredFields = {'SD', 't', 'd'};
    for i = 1:length(requiredFields)
        assert(isfield(data, requiredFields{i}), ['Missing required field: ', requiredFields{i}]);
    end
    assert(isfield(data.SD, 'Lambda'), 'Missing required field: Lambda in SD');
    assert(isnumeric(data.SD.Lambda) && numel(data.SD.Lambda) == 2, 'Lambda should be a numeric array with 2 elements');
    
    extCoeffData = readtable(extinctionCoefficientsFile);
    DPFData = readtable(DPFperTissueFile, 'Format', '%s%f');
    relDPFData = readtable(relDPFfile);

    Lambda = data.SD.Lambda;
    time = data.t;
    intensityData = data.d;
    nChannels = 20;

    % Validating size and type of intensityData
    assert(isnumeric(intensityData) && size(intensityData, 2) == 2 * nChannels, 'Intensity data should be a numeric matrix with 2*nChannels columns');

    % Validating tissueType exists in DPFData
    tissueIdx = find(strcmp(DPFData{:, 1}, tissueType));
    assert(~isempty(tissueIdx), 'Tissue type not found in DPF per tissue file');
    
    DPF = DPFData{tissueIdx, 2};
    relDPF1 = interp1(relDPFData{:, 1}, relDPFData{:, 2}, Lambda(1));
    relDPF2 = interp1(relDPFData{:, 1}, relDPFData{:, 2}, Lambda(2));

    % Validating wavelengths exist in extCoeffData
    assert(any(extCoeffData{:, 1} == Lambda(1)), 'Lambda(1) not found in extinction coefficients file');
    assert(any(extCoeffData{:, 1} == Lambda(2)), 'Lambda(2) not found in extinction coefficients file');

    % Calculating the effective L
    Leff1 = DPF * SDS * relDPF1;
    Leff2 = DPF * SDS * relDPF2;

    % Extinction coefficients for the two wavelengths
    extCoeff1 = extCoeffData{extCoeffData{:, 1} == Lambda(1), 3:4};
    extCoeff2 = extCoeffData{extCoeffData{:, 1} == Lambda(2), 3:4};

    % Calculating ODs
    OD1 = -log10(intensityData(:, 1:nChannels) ./ intensityData(1, 1:nChannels));
    OD2 = -log10(intensityData(:, nChannels+1:end) ./ intensityData(1, nChannels+1:end));

    % Calculating ΔHbR and ΔHbO
    invExtCoeff = inv([extCoeff1; extCoeff2]);
    dHb = zeros(length(time), nChannels, 2);

    for ch = 1:nChannels
        for t = 1:length(time)
            dHb(t, ch, :) = invExtCoeff * [OD1(t, ch) / Leff1; OD2(t, ch) / Leff2];
        end
    end

    dHbO = squeeze(dHb(:, :, 1));
    dHbR = squeeze(dHb(:, :, 2));

    % Plotting ΔHbR and ΔHbO for specified channels
    fig = [];
    if ~isempty(plotChannelIdx)
        fig = figure;
        for ch = plotChannelIdx
            subplot(length(plotChannelIdx), 1, find(plotChannelIdx == ch));
            plot(time, dHbO(:, ch), 'r', time, dHbR(:, ch), 'b');
            title(['Channel ' num2str(ch)]);
            xlabel('Time (s)');
            ylabel('Concentration Change (μM)');
            legend('ΔHbO', 'ΔHbR');
        end
    end

    % Computing FFT and SNR only for the first channel from the first file
    if strcmp(dataFile, 'FN_031_V2_Postdose2_Nback.mat')
        fs = 1 / mean(diff(time)); % Sampling frequency
        n = length(time);
        
        % Computing the FFTs 
        Y_O = fft(dHbO(:, 1));
        P2_O = abs(Y_O / n);
        P1_O = P2_O(1:floor(n/2)+1);
        P1_O(2:end-1) = 2*P1_O(2:end-1);
        
        Y_R = fft(dHbR(:, 1));
        P2_R = abs(Y_R / n);
        P1_R = P2_R(1:floor(n/2)+1);
        P1_R(2:end-1) = 2*P1_R(2:end-1);
        
        f = fs * (0:(n/2)) / n;

        % Heartbeat frequency peak
        [~, peakIdx] = max(P1_O(f >= 1 & f <= 2));
        heartbeatFreq = f(peakIdx);

        % Calculating the SNRs
        signal_O = P1_O(peakIdx);
        noise_O = mean(P1_O(f > 2.5));
        signal_R = P1_R(peakIdx);
        noise_R = mean(P1_R(f > 2.5));

        SNR_O = signal_O / noise_O;
        disp(['First SNR (ΔHbO) is: ', num2str(SNR_O)]);
        SNR_R = signal_R / noise_R;
        disp(['Second SNR (ΔHbR) is: ', num2str(SNR_R)]);

        % Plotting the Fourier transform 
        figure;
        plot(f, P1_O, 'r', f, P1_R, 'b');
        title('Fourier transform of the first channel [ΔHbO and ΔHbR]');
        xlabel('Frequency (Hz)');
        ylabel('|P1(f)|');
        legend('ΔHbO', 'ΔHbR');
    end
end