% Import CSV File
clearvars;
% Write HyperParameter
File_PATH = 'Data/Test';
Result_Path='Data/Result';
Select_variable='Test';

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

% Concat The data and reshpe [3 * 174 * length]

ABS_Data=cat(3,abs(cfs_Id),abs(cfs_Vds),abs(cfs_Vgs));
ABS_Data= reshape(ABS_Data, 3,[]);


angle_Data=cat(3,angle(cfs_Id),angle(cfs_Vds),angle(cfs_Vgs));
angle_Data= reshape(angle_Data, 3,[]);


% Save csv file it can't save multi dim file in matlab code so it need to
% unsquuez file when file import
test=abs(cfs_Id(1:3,1:10));
test=reshape(test,3,2,5);
writematrix(test,strcat(Result_Path,'/Wavelet_',Select_variable,'abs.csv'));
%writematrix(angle_Data,strcat(Result_Path,'/Wavelet_',Select_variable,'angle.csv'));