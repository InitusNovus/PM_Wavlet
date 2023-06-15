clearvars;
%% Global Constant Definitions
ConstGlb = struct;
ConstGlb.whatTrs = ["turnon", "turnoff"];
ConstGlb.whatCh = ["Id", "Vds", "Vgs"];
ConstGlb.testInductance = 1.0;  % [mH]
ConstGlb.namePart = dictionary( ...
    ["Ron" "Roff" "Pulse" "Vds" "Vgson" "Vgsoff" "trs" "id"], ...
    [ 2     4      6       8     10      12       13    16]);
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
            data.(whatTrs)(k).(whatCh) = nan(L_h, L_w);
            data.(whatTrs)(k).(whatCh) = dataRaw.(whatTrs)(L_h*(i-1) + 1: L_h*(i), :, k);
            
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























