clearvars;
%% Global Constant Definitions
ConstGlb = struct;
ConstGlb.whatTrs = ["turnon", "turnoff"];
ConstGlb.whatCh = ["Id", "Vds", "Vgs"];
ConstGlb.testInductance = 1.0;  % [mH]
ConstGlb.namePart = dictionary( ...
    ["Ron" "Roff" "Pulse" "Vds" "Vgson" "Vgsoff" "trs" "id"], ...
    [ 2     4      6       8     10      12       13    16]);
ConstGlb.minmax = dictionary( ...
    ["Ronmin"   "Ronmax"    "Roffmin"   "Roffmax"   "Idmin"     "Idmax" ...
     "Vdsmin"   "Vdsmax"    "Vgsonmin"  "Vgsonmax"  "Vgsoffmin" "Vgsoffmax"], ...
    [0          30          0           30          0           36 ...
     200        600         10          20          -5          0]);
ConstGlb.dataSize = [87 250];

%% Import Data
dataImportTimer = tic;
dataPath = "../Data/Wavelet_Result_Transient/Wavelet";
dataRaw = struct;
for whatTrs = ConstGlb.whatTrs
    L_dataPath = dataPath + "/" + whatTrs;
    L_fileList = dir(L_dataPath);
    L_fileList = L_fileList(3:end); % Remove the first two files: "." and ".."
    L_numFile = numel(L_fileList);
    dataRaw.(whatTrs) = nan([ConstGlb.dataSize L_numFile] .* [6 1 1]);
    for i = 1:L_numFile
        L_fileFullName = L_dataPath + "/" + L_fileList(i).name;
        dataRaw.(whatTrs)(:,:,i) = readmatrix(L_fileFullName);
    end
    clear i L_fileFullName
    dataList.(whatTrs) = L_fileList;
end
clear whatTrs L_dataPath L_fileList L_numFile
timeDataImport = toc(dataImportTimer);
fprintf("Data Import Time: %02d:%02d:%02.3f\n", sec2time(timeDataImport));
clear dataImportTimer

%% Parse the Data
dataParseTimer = tic;
data = struct;
L_whatTrs = string(fieldnames(dataRaw));
for j = 1:length(L_whatTrs)
    whatTrs = L_whatTrs(j);
    L_numData = size(dataRaw.(whatTrs), 3);
    L_h = ConstGlb.dataSize(1);
    L_w = ConstGlb.dataSize(2);
    for k = 1:L_numData
        for i = 1:length(ConstGlb.whatCh)
            whatCh = ConstGlb.whatCh(i);
            data.(whatTrs)(k).(whatCh).abs = nan(L_h, L_w);
            data.(whatTrs)(k).(whatCh).abs = dataRaw.(whatTrs)(L_h*(i-1) + 1: L_h*(i), :, k);
            data.(whatTrs)(k).(whatCh).angle = nan(L_h, L_w);
            data.(whatTrs)(k).(whatCh).angle = dataRaw.(whatTrs)(L_h*(i-1 + 3) + 1: L_h*(i + 3), :, k);
            
            % Parse the labels
            L_fileName = dataList.(whatTrs)(k).name;
            L_nameParts = split(L_fileName, "_");
            dict = ConstGlb.namePart;
            Ron = string(L_nameParts{dict("Ron")});
            Ron = double(Ron)/10;      % [Ohm]
            Roff = string(L_nameParts{dict("Roff")});
            Roff = double(Roff)/10;    % [Ohm]
            Pulse = string(L_nameParts{dict("Pulse")});
            Pulse = double(Pulse);     % [us]
            Vds = string(L_nameParts{dict("Vds")});
            Vds = double(Vds);         % [V]
            Vgson = string(L_nameParts{dict("Vgson")});
            Vgson = double(Vgson);     % [V]
            Vgsoff = string(L_nameParts{dict("Vgsoff")});
            Vgsoff = -double(Vgsoff);   % [V]
            Id = Vds/(ConstGlb.testInductance*1e-3) * (Pulse*1e-6);
            data.(whatTrs)(k).labels = dictionary( ...
                ["Ron" "Roff" "Id" "Vds" "Vgson" "Vgsoff"], ...
                [ Ron   Roff   Id   Vds   Vgson   Vgsoff]);
        end
        clear i whatCh
        clear dict Ron Roff Pulse Vds Vgson Vgsoff Id
        clear L_fileName L_nameParts
    end
    clear k
end
clear j L_whatTrs whatTrs L_numData L_h L_w
timeDataParse = toc(dataParseTimer);
fprintf("Data Parse Time: %02d:%02d:%02.3f\n", sec2time(timeDataParse));
clear dataParseTimer

%% Parse the Train and Test Data
whatTrs = "turnon";
whatCh = "Vds";
L_data = data.(whatTrs);
numData = length(data.(whatTrs));
numInputs = length(keys(L_data(1).labels));

randomSeed = 917;

% random seed
rng(randomSeed);

% predictors
L_data = data.(whatTrs);
L_minmax = ConstGlb.minmax;
predictors = array2table(nan(numData, numInputs));
predictors.Properties.VariableNames = keys(L_data(1).labels);

responses = nan([ConstGlb.dataSize 2 numData]);
for k = 1:numData
    predictors{k,:} = values(L_data(k).labels)';
    responses(:,:,1,k) = L_data(k).(whatCh).abs;
    responses(:,:,2,k) = L_data(k).(whatCh).angle;
end

% Normalization
for vars = string(predictors.Properties.VariableNames)
    predictors.(vars) = ...
        (predictors.(vars) - L_minmax(vars+"min"))./(L_minmax(vars+"max") - L_minmax(vars+"min"));
end
responses(:,:,1,:) = responses(:,:,1,:)./10;
responses(:,:,2,:) = (responses(:,:,2,:) + pi)./(2*pi);

% Generate a permutation of indices
shuffleIdx = randperm(numData);

% Shuffle predictors and responses
predictors = predictors(shuffleIdx, :);
responses = responses(:, :, :, shuffleIdx);

% Define the proportion of data to be used for training
trainProportion = 0.8;  % or whatever proportion you prefer

% Calculate the number of training samples
numTrain = round(trainProportion * numData);

% Training set
XTrain = predictors{1:numTrain, :};
XTrain = reshape(XTrain', 1, 1, numInputs, length(XTrain));
YTrain = responses(:, :, :, 1:numTrain);

XTrain = single(XTrain);
YTrain = single(YTrain);

% Validation set
XVal = predictors{numTrain+1:end, :};
XVal = reshape(XVal', 1, 1, numInputs, length(XVal));
YVal = responses(:, :, :, numTrain+1:end);

XVal = single(XVal);
YVal = single(YVal);

%% Generator

dropoutTConv = 0.5;
dropoutFcl = 0.5;
dropoutResize = 0.5;

% The definition of the layer graph
layers1 = [
    % featureInputLayer(6, 'Normalization', 'none', 'Name', 'input') % assuming your input vector has 6 elements
    imageInputLayer([1 1 numInputs], 'Normalization','none','Name', 'input')

    fullyConnectedLayer(256)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    fullyConnectedLayer(512)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl, 'Name', 'dropout fcl')];

layers2_1 = [ 
    fullyConnectedLayer(512, 'Name', 'fcl 1_1')
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    fullyConnectedLayer(512)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    resize3dLayer("OutputSize", [1 64 8], 'Method','trilinear', 'Name', 'resize 1_1')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 1_1')
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 8)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [3 192 4], 'Method','trilinear', 'Name', 'resize 2_1')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 2_1')
    dropoutLayer(dropoutResize, 'Name', 'dropout resize 2_1')];

layers2_2 = [ 
    fullyConnectedLayer(512, 'Name', 'fcl 1_2')
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    fullyConnectedLayer(512)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    resize3dLayer("OutputSize", [1 64 8], 'Method','trilinear', 'Name', 'resize 1_2')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 1_2')
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 8)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [3 192 4], 'Method','trilinear', 'Name', 'resize 2_2')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 2_2')
    dropoutLayer(dropoutResize, 'Name', 'dropout resize 2_2')];

layers2_3 = [ 
    fullyConnectedLayer(512, 'Name', 'fcl 1_3')
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    fullyConnectedLayer(512)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)
    
    resize3dLayer("OutputSize", [1 64 8], 'Method','trilinear', 'Name', 'resize 1_3')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 1_3')
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 8)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [3 192 4], 'Method','trilinear', 'Name', 'resize 2_3')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 2_3')
    dropoutLayer(dropoutResize, 'Name', 'dropout resize 2_3')];

layers2_4 = [ 
    fullyConnectedLayer(512, 'Name', 'fcl 1_4')
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    fullyConnectedLayer(512)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)
    
    resize3dLayer("OutputSize", [1 64 8], 'Method','trilinear', 'Name', 'resize 1_4')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 1_4')
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 8)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [3 192 4], 'Method','trilinear', 'Name', 'resize 2_4')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 2_4')
    dropoutLayer(dropoutResize, 'Name', 'dropout resize 2_4')];

layers2_5 = [ 
    fullyConnectedLayer(512, 'Name', 'fcl 1_5')
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    fullyConnectedLayer(512)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)
    
    resize3dLayer("OutputSize", [1 64 8], 'Method','trilinear', 'Name', 'resize 1_5')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 1_5')
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 8)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [3 192 4], 'Method','trilinear', 'Name', 'resize 2_5')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 2_5')
    dropoutLayer(dropoutResize, 'Name', 'dropout resize 2_5')];

layers2_6 = [ 
    fullyConnectedLayer(512, 'Name', 'fcl 1_6')
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    fullyConnectedLayer(512)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)
    
    resize3dLayer("OutputSize", [1 64 8], 'Method','trilinear', 'Name', 'resize 1_6')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 1_6')
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 8)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [3 192 4], 'Method','trilinear', 'Name', 'resize 2_6')
    batchNormalizationLayer();
    leakyReluLayer('Name', 'activation resize 2_6')
    dropoutLayer(dropoutResize, 'Name', 'dropout resize 2_6')];

layers3 = [
    concatenationLayer(1,6,'Name','concat')
    
    transposedConv2dLayer([7 7], 128)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([7 7], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([7 7], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([7 7], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [160 160 4], 'Method','trilinear', 'Name', 'resize 3')
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 63)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [75 238 8], 'Method','trilinear', 'Name', 'resize 4')
    batchNormalizationLayer();
    leakyReluLayer()

    transposedConv2dLayer([5 5], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 16)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [87 250 2], 'Method','trilinear', 'Name', 'resize final')
    
    regressionLayer('Name','regressionoutput')];

lgraph = layerGraph(layers1);
lgraph = addLayers(lgraph, layers2_1);
lgraph = connectLayers(lgraph, 'dropout fcl','fcl 1_1');
lgraph = addLayers(lgraph, layers2_2);
lgraph = connectLayers(lgraph, 'dropout fcl', 'fcl 1_2');
lgraph = addLayers(lgraph, layers2_3);
lgraph = connectLayers(lgraph, 'dropout fcl', 'fcl 1_3');
lgraph = addLayers(lgraph, layers2_4);
lgraph = connectLayers(lgraph, 'dropout fcl', 'fcl 1_4');
lgraph = addLayers(lgraph, layers2_5);
lgraph = connectLayers(lgraph, 'dropout fcl', 'fcl 1_5');
lgraph = addLayers(lgraph, layers2_6);
lgraph = connectLayers(lgraph, 'dropout fcl', 'fcl 1_6');
lgraph = addLayers(lgraph, layers3);
lgraph = connectLayers(lgraph, 'dropout resize 2_1', 'concat/in1');
lgraph = connectLayers(lgraph, 'dropout resize 2_2', 'concat/in2');
lgraph = connectLayers(lgraph, 'dropout resize 2_3', 'concat/in3');
lgraph = connectLayers(lgraph, 'dropout resize 2_4', 'concat/in4');
lgraph = connectLayers(lgraph, 'dropout resize 2_5', 'concat/in5');
lgraph = connectLayers(lgraph, 'dropout resize 2_6', 'concat/in6');

% plot(lgraph)

% lgraph = layerGraph(layers);
analyzeNetwork(lgraph)
%% Training 
% Training options
maxEpochs = 50000;
miniBatchSize = 32;
validationFreq = 10;

options = trainingOptions('sgdm', ...
    'MaxEpochs',maxEpochs, ...
    'L2Regularization', 0.0001, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',1000, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress', ...
    'Shuffle','every-epoch', ...
    'MiniBatchSize', miniBatchSize, ...
    'ExecutionEnvironment', 'gpu', ...
    'ValidationData', {XVal, YVal}, ...
    'ValidationFrequency', validationFreq, ...
    'OutputNetwork', 'last-iteration');

% Train the network
[net, info] = trainNetwork(XTrain, YTrain, lgraph, options);

%% 
k = randi(size(XVal, 4));
YPred = predict(net, XVal(1,1,:,k));

figure(1);
imagesc(YPred(:,:,1)); % Plot the CWT coefficients using imagesc
colormap(jet); % Choose the colormap for the CWT coefficients plot
colorbar; % Show the color bar for the CWT coefficients
xlabel('Time');
ylabel('Frequency');
% set(gca, 'YScale', 'log'); % Set the y-axis to log scale to better visualize the frequency components
axis('tight')
title("Wavelet Data 250 x 87 (Abs), Predict")

figure(2);
imagesc(YVal(:,:,1,k)); % Plot the CWT coefficients using imagesc
colormap(jet); % Choose the colormap for the CWT coefficients plot
colorbar; % Show the color bar for the CWT coefficients
xlabel('Time');
ylabel('Frequency');
% set(gca, 'YScale', 'log'); % Set the y-axis to log scale to better visualize the frequency components
axis('tight')
title("Wavelet Data 250 x 87 (Abs), GT")

figure(3);
imagesc(YPred(:,:,2)); % Plot the CWT coefficients using imagesc
colormap(jet); % Choose the colormap for the CWT coefficients plot
colorbar; % Show the color bar for the CWT coefficients
xlabel('Time');
ylabel('Frequency');
% set(gca, 'YScale', 'log'); % Set the y-axis to log scale to better visualize the frequency components
axis('tight')
title("Wavelet Data 250 x 87 (Angle), Predict")

figure(4);
imagesc(YVal(:,:,2,k)); % Plot the CWT coefficients using imagesc
colormap(jet); % Choose the colormap for the CWT coefficients plot
colorbar; % Show the color bar for the CWT coefficients
xlabel('Time');
ylabel('Frequency');
% set(gca, 'YScale', 'log'); % Set the y-axis to log scale to better visualize the frequency components
axis('tight')
title("Wavelet Data 250 x 87 (Angle), GT")

%
reconCfs.abs = double(YPred(:,:,1));
reconCfs.angle = double(YPred(:,:,2));
maxRawAngle = max(max(reconCfs.angle));
minRawAngle = min(min(reconCfs.angle));
reconCfs.angle = (reconCfs.angle - minRawAngle) * (2*pi/(maxRawAngle - minRawAngle)) - pi;
reconCfs.complex = reconCfs.abs .* exp(complex(0,1) * reconCfs.angle);
reconCfs.waveform = icwt(reconCfs.complex);

figure(5); 
plot(reconCfs.waveform)
grid on;















