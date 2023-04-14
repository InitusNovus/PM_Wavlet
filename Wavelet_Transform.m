%% Import CSV File
clearvars;
% Write HyperParameter
File_PATH = 'Data/Test';
Result_Path='Data/Result';
Select_variable='Test';

File_List=dir(File_PATH);
[Data_num,~] = size(File_List);

VOP=5;

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
[cfs_Id, f_Id] = cwt(Id, 'Morse', fs ,VoicePerOctave=VOP);
[cfs_Vds, f_Vds] = cwt(Vds, 'Morse', fs ,VoicePerOctave=VOP);
[cfs_Vgs, f_Vgs] = cwt(Vgs, 'Morse', fs ,VoicePerOctave=VOP);

% Resize Data

original_num_columns = len(cfs_Id);

num_columns_resampled = len(cfs_Id)/50; % New number of columns
resample_factor = size(cfs_Id, 2) / num_columns_resampled;
t_res = resample(t.', 1, resample_factor).';

% ID
cfs_Id_res = resample(cfs_Id.', 1, resample_factor).';

% VDS
cfs_Vds_res = resample(cfs_Vds.', 1, resample_factor).';

% VGS
cfs_Vgs_res = resample(cfs_Vgs.', 1, resample_factor).';



% % Concat The data and reshpe [3 * 174 * length]
% 
% ABS_Data=cat(3,abs(cfs_Id),abs(cfs_Vds),abs(cfs_Vgs));
% ABS_Data= reshape(ABS_Data, 3,[]);
% 
% 
% angle_Data=cat(3,angle(cfs_Id),angle(cfs_Vds),angle(cfs_Vgs));
% angle_Data= reshape(angle_Data, 3,[]);
% 
% 
% % Save csv file it can't save multi dim file in matlab code so it need to
% % unsquuez file when file import
% test=abs(cfs_Id(1:3,1:10));
% test=reshape(test,3,2,5);
% writematrix(test,strcat(Result_Path,'/Wavelet_',Select_variable,'abs.csv'));
% %writematrix(angle_Data,strcat(Result_Path,'/Wavelet_',Select_variable,'angle.csv'));