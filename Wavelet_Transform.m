% Import CSV File
clearvars;
% Write HyperParameter
File_PATH = 'Data/Test';

File_List=dir(File_PATH);
[Data_num,~] = size(File_List);

% File Selection Section

file_name=File_List(3).name;
file_name=strcat(File_PATH,'/',file_name);
data = readtable(file_name);

% Data parsing

data = renamevars(data, "TIME", "Time");
t = data.Time;
Vgs = data.Vgs;
Vds = data.Vds;
Id = data.Id;

% Perform Contimuous Wavelet Transform

fs = length(t)/(t(length(t)) - t(1));
[cfs_Id, f_Id] = cwt(Id, 'Morse', fs);
[cfs_Vds, f_Vds] = cwt(Vds, 'Morse', fs);
[cfs_Vgs, f_Vgs] = cwt(Vgs, 'Morse', fs);




