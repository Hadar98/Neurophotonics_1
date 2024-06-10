%% Hadar Granot 208749176 and Maya Wertsman 208915355

% Inputs
dataFile1 = 'FN_031_V2_Postdose2_Nback.mat';
dataFile2 = 'FN_032_V1_Postdose1_Nback.mat';
SDS = 3; % in cm
tissueType = 'adult_head'; % For example
plotChannelIdx = [1, 2];
extCoeffFile = 'ExtinctionCoefficientsData.csv';
DPFFile = 'DPFperTissue.txt';
relDPFfile = 'RelativeDPFCoefficients.csv';

% Running the function for the first file
[dHbR1, dHbO1, fig1] = CalcNIRS(dataFile1, SDS, tissueType, plotChannelIdx, extCoeffFile, DPFFile, relDPFfile);

% Running the function for the second file
[dHbR2, dHbO2, fig2] = CalcNIRS(dataFile2, SDS, tissueType, plotChannelIdx, extCoeffFile, DPFFile, relDPFfile);
