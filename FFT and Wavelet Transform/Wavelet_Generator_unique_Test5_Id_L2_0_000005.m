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
whatCh = "Id";
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

% Group data by unique labels
[uniqueLabels, ~, ic] = unique(predictors, 'rows');
groupedIndices = accumarray(ic, 1:size(predictors, 1), [], @(x) {x});

selectedData = [];

for idx = 1:length(groupedIndices)
    group = groupedIndices{idx};
    
    % Randomly shuffle the group
    numInGroup = length(group);
    shuffledGroup = group(randperm(numInGroup));
    
    selectedData = [selectedData; shuffledGroup(1)];
    % Training set for this group
    % trainIndices = shuffledGroup(1:numTrainInGroup);
    % XTrainGroup = predictors{trainIndices, :};
    % XTrainGroup = reshape(XTrainGroup', 1, 1, numInputs, length(XTrainGroup));
    % YTrainGroup = responses(:, :, :, trainIndices);
    % 
    % XTrain = cat(4, XTrain, single(XTrainGroup));
    % YTrain = cat(4, YTrain, single(YTrainGroup));
    % 
    % % Validation set for this group
    % if numTrainInGroup < numInGroup  % Only if there's data left for validation
    %     valIndices = shuffledGroup(numTrainInGroup+1:end);
    %     XValGroup = predictors{valIndices, :};
    %     XValGroup = reshape(XValGroup', 1, 1, numInputs, length(XValGroup));
    %     YValGroup = responses(:, :, :, valIndices);
    % 
    %     XVal = cat(4, XVal, single(XValGroup));
    %     YVal = cat(4, YVal, single(YValGroup));
    % end
end

predictors = predictors(selectedData, :);
responses = responses(:, :, :, selectedData);
numData = size(predictors, 1);

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

dropoutTConv = 0.3;
dropoutFcl = 0.5;
dropoutResize = 0.5;

% The definition of the layer graph
layers1 = [
    % featureInputLayer(6, 'Normalization', 'none', 'Name', 'input') % assuming your input vector has 6 elements
    imageInputLayer([1 1 numInputs], 'Normalization','none','Name', 'input')

    % fullyConnectedLayer(512)
    % batchNormalizationLayer();
    % leakyReluLayer()
    % % tanhLayer()
    % dropoutLayer(dropoutFcl)

    fullyConnectedLayer(30000)
    batchNormalizationLayer();
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl, 'Name', 'dropout fcl')];

layers2_1 = createTConvTall(512, dropoutFcl, dropoutTConv, dropoutResize, "layers2_1");
layers2_2 = createTConvTall(512, dropoutFcl, dropoutTConv, dropoutResize, "layers2_2");
layers2_3 = createTConvTall(512, dropoutFcl, dropoutTConv, dropoutResize, "layers2_3");
layers2_4 = createTConvTall(512, dropoutFcl, dropoutTConv, dropoutResize, "layers2_4");
layers2_5 = createTConvTall(512, dropoutFcl, dropoutTConv, dropoutResize, "layers2_5");
layers2_6 = createTConvTall(512, dropoutFcl, dropoutTConv, dropoutResize, "layers2_6");

layersConcatTall = [
    concatenationLayer(1,6,'Name','concatTall')];

layers3_1 = createTConvFat(512, dropoutFcl, dropoutTConv, dropoutResize, "layers3_1");
layers3_2 = createTConvFat(512, dropoutFcl, dropoutTConv, dropoutResize, "layers3_2");
layers3_3 = createTConvFat(512, dropoutFcl, dropoutTConv, dropoutResize, "layers3_3");
layers3_4 = createTConvFat(512, dropoutFcl, dropoutTConv, dropoutResize, "layers3_4");
layers3_5 = createTConvFat(512, dropoutFcl, dropoutTConv, dropoutResize, "layers3_5");
layers3_6 = createTConvFat(512, dropoutFcl, dropoutTConv, dropoutResize, "layers3_6");

layersConcatFat = [
    concatenationLayer(2,6,'Name','concatFat')];

layersConcatOne = [
    depthConcatenationLayer(2, "Name", "concatOne")];

layers4_1 = createTConvFinal(dropoutFcl, dropoutTConv, dropoutResize, "layers4_1");
layers4_2 = createTConvFinal(dropoutFcl, dropoutTConv, dropoutResize, "layers4_2");

layersConcatFinal = [
    depthConcatenationLayer(2, "Name", "concatFinal")];

layersOut = [
    regressionLayer('Name','regressionoutput')];

lgraph = layerGraph(layers1);
lgraph = addLayers(lgraph, layers2_1);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers2_1_in');
lgraph = addLayers(lgraph, layers2_2);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers2_2_in');
lgraph = addLayers(lgraph, layers2_3);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers2_3_in');
lgraph = addLayers(lgraph, layers2_4);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers2_4_in');
lgraph = addLayers(lgraph, layers2_5);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers2_5_in');
lgraph = addLayers(lgraph, layers2_6);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers2_6_in');
lgraph = addLayers(lgraph, layersConcatTall);
lgraph = connectLayers(lgraph, 'layers2_1_out', 'concatTall/in1');
lgraph = connectLayers(lgraph, 'layers2_2_out', 'concatTall/in2');
lgraph = connectLayers(lgraph, 'layers2_3_out', 'concatTall/in3');
lgraph = connectLayers(lgraph, 'layers2_4_out', 'concatTall/in4');
lgraph = connectLayers(lgraph, 'layers2_5_out', 'concatTall/in5');
lgraph = connectLayers(lgraph, 'layers2_6_out', 'concatTall/in6');

lgraph = addLayers(lgraph, layers3_1);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers3_1_in');
lgraph = addLayers(lgraph, layers3_2);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers3_2_in');
lgraph = addLayers(lgraph, layers3_3);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers3_3_in');
lgraph = addLayers(lgraph, layers3_4);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers3_4_in');
lgraph = addLayers(lgraph, layers3_5);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers3_5_in');
lgraph = addLayers(lgraph, layers3_6);
lgraph = connectLayers(lgraph, 'dropout fcl', 'layers3_6_in');
lgraph = addLayers(lgraph, layersConcatFat);
lgraph = connectLayers(lgraph, 'layers3_1_out', 'concatFat/in1');
lgraph = connectLayers(lgraph, 'layers3_2_out', 'concatFat/in2');
lgraph = connectLayers(lgraph, 'layers3_3_out', 'concatFat/in3');
lgraph = connectLayers(lgraph, 'layers3_4_out', 'concatFat/in4');
lgraph = connectLayers(lgraph, 'layers3_5_out', 'concatFat/in5');
lgraph = connectLayers(lgraph, 'layers3_6_out', 'concatFat/in6');

lgraph = addLayers(lgraph, layersConcatOne);
lgraph = connectLayers(lgraph, 'concatTall/out', 'concatOne/in1');
lgraph = connectLayers(lgraph, 'concatFat/out', 'concatOne/in2');

lgraph = addLayers(lgraph, layers4_1);
lgraph = connectLayers(lgraph, 'concatOne/out', 'layers4_1_in');
lgraph = addLayers(lgraph, layers4_2);
lgraph = connectLayers(lgraph, 'concatOne/out', 'layers4_2_in');

lgraph = addLayers(lgraph, layersConcatFinal);
lgraph = connectLayers(lgraph, 'layers4_1_out', 'concatFinal/in1');
lgraph = connectLayers(lgraph, 'layers4_2_out', 'concatFinal/in2');

lgraph = addLayers(lgraph, layersOut);
lgraph = connectLayers(lgraph, 'concatFinal/out', 'regressionoutput');


% plot(lgraph)

% lgraph = layerGraph(layers);
analyzeNetwork(lgraph)
%% Training
% Training options
% maxEpochs = 50000;
maxEpochs = 5000;
miniBatchSize = 16;
validationFreq = 22;
L2Reg = 0.000005;

% options = trainingOptions('sgdm', ...
%     'MaxEpochs',maxEpochs, ...
%     'L2Regularization', 0.0001, ...
%     'GradientThreshold',1, ...
%     'InitialLearnRate',0.01, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropPeriod',1000, ...
%     'LearnRateDropFactor',0.1, ...
%     'Verbose',0, ...
%     'Plots','training-progress', ...
%     'Shuffle','every-epoch', ...
%     'MiniBatchSize', miniBatchSize, ...
%     'ExecutionEnvironment', 'gpu', ...
%     'ValidationData', {XVal, YVal}, ...
%     'ValidationFrequency', validationFreq, ...
%     'OutputNetwork', 'last-iteration');

options = trainingOptions('sgdm', ...
    'MaxEpochs',maxEpochs, ...
    'L2Regularization', L2Reg, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.05, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',1000, ...
    'LearnRateDropFactor',0.5, ...
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

referenceNum = shuffleIdx(numTrain+1);
referencePath = "../Data/Wavelet_Result_Transient/TimeDomain" + "/" + whatTrs;
referenceList = dir(referencePath);
referenceList = referenceList(3:end);
referenceFileName = (referencePath + "/" + string(referenceList(referenceNum).name));
referenceData = readtable(referenceFileName); 

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

plotRecon = normalize(reconCfs.waveform);
plotT = (1:250)*(12500/250);
plotRef = normalize(referenceData.(whatCh));

figure(5);
% plot(reconCfs.waveform)
plot(plotRef)
hold on;
plot(plotT, plotRecon, LineWidth=1.5)
hold off;
grid on;
axis("tight")
legend("Actual", "Predicted", "Location", "best")
title("Waveform Reconstructed from Predicted Wavelet Image", FontWeight="bold")
xlabel("Time step (of reference waveform)")
ylabel("Amplitude (normalized)")


%% Function Definitions

function layers = createTConvTall(widthFcl, dropoutFcl, dropoutTConv, dropoutResize, namePrefix)

namePrefix = string(namePrefix);
layers = [
    fullyConnectedLayer(widthFcl, "Name", namePrefix + "_in")
    batchNormalizationLayer()
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    % fullyConnectedLayer(widthFcl)
    % batchNormalizationLayer()
    % leakyReluLayer()
    % % tanhLayer()
    % dropoutLayer(dropoutFcl)

    resize3dLayer("OutputSize", [4 64 128], 'Method','trilinear')
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [32 192 24], "Method","trilinear")
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutResize, "Name", namePrefix + "_out")];

end

function layers = createTConvFat(widthFcl, dropoutFcl, dropoutTConv, dropoutResize, namePrefix)

namePrefix = string(namePrefix);
layers = [
    fullyConnectedLayer(widthFcl, "Name", namePrefix + "_in")
    batchNormalizationLayer()
    leakyReluLayer()
    % tanhLayer()
    dropoutLayer(dropoutFcl)

    % fullyConnectedLayer(widthFcl)
    % batchNormalizationLayer()
    % leakyReluLayer()
    % % tanhLayer()
    % dropoutLayer(dropoutFcl)

    resize3dLayer("OutputSize", [64 4 128], 'Method','trilinear')
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer()
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 128)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [192 32 24], "Method","trilinear")
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutResize, "Name", namePrefix + "_out")];

end


function layers = createTConvFinal(dropoutFcl, dropoutTConv, dropoutResize,namePrefix)
layers = [
    transposedConv2dLayer([7 7], 64, "Name", namePrefix + "_in")
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([7 7], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([7 7], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([7 7], 64)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [160 160 16], 'Method','trilinear')
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutResize)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    transposedConv2dLayer([5 5], 32)
    batchNormalizationLayer();
    leakyReluLayer()
    dropoutLayer(dropoutTConv)

    resize3dLayer("OutputSize", [75 238 16], 'Method','trilinear')
    batchNormalizationLayer();
    leakyReluLayer()

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

    resize3dLayer("OutputSize", [87 250 1], 'Method','trilinear', 'Name', namePrefix + "_out")];
end


















