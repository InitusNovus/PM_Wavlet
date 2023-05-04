%% TODO:
%   Make double -> single to reduce the size of result data
%   The precision is not important.
%% Import CSV File
clearvars;
% Write HyperParameter
File_PATH = 'Data/Ron_50_Roff_50';
Result_Path='Data/FFT_toSingle/';
Select_variable='Test';
% Please write Resist
R_on=50;
R_off=50;

File_List=dir(File_PATH);
[Data_num,~] = size(File_List);

VPO=5;
%% File Selection Section
totalTimer = tic;
for i = 3: Data_num
    iterTimer = tic;

    file_name=File_List(i).name;
    name_parts=split(file_name,'_');
    file_name=strcat(File_PATH,'/',file_name);
    data = readtable(file_name);
    % Data parsing
    
    data = renamevars(data, "Labels", "Time");
    t = data.Time;
    Vgs = data.Vgs;
    Vds = data.Vds;
    Id = data.Id;
    
    % Eliminate the mean value in Vds to finde the switching transient 

    Vds_mean = mean(Vds);
    Vds_zeromean = Vds - Vds_mean;
    min_gap = 10000; % Arbitrary value

    zero_crossings = Find_zeroCrossing(Vds_zeromean, min_gap);
    assert(length(zero_crossings) == 4, "Error: Switching transient detection failed")

    % Interpolate the time values corresponding to the zero crossings
    zero_crossing_times = interp1(1:length(t), t, zero_crossings);
    clear zero_crossings

    zero_crossing_names = {"Charging"; "Turn-Off"; "Turn-On"; "End"}
    st_time = containers.Map();
    for i = 1:4
    st_time(zero_crossing_names{i}) = zero_crossing_times(i);
    end
    clear i

    st_time_range = 0.1 * [1e-6 11e-6];
    
    %Turn-OFF transient
    turnoff_time_range = st_time("Turn-Off") + st_time_range;
    index_turnoff = (turnoff_time_range(1) < t) & (t < turnoff_time_range(2));
    t_turnoff = t(index_turnoff);
    Vgs_turnoff = Vgs(index_turnoff);
    Vds_turnoff = Vds(index_turnoff);
    Id_turnoff = Id(index_turnoff);

    %Turn-On Transient
    turnon_time_range = st_time("Turn-On") + st_time_range;
    index_turnon = (turnon_time_range(1) < t) & (t < turnon_time_range(2));
    t_turnon = t(index_turnon);
    Vgs_turnon = Vgs(index_turnon);
    Vds_turnon = Vds(index_turnon);
    Id_turnon = Id(index_turnon);

    % clear

    clear turnoff_time_range turnon_time_range
    clear index_turnoff index_turnon
    clear st_time_range

    % Data to table

    data_turnoff = table();
    data_turnoff.t = t_turnoff;
    data_turnoff.Vgs = Vgs_turnoff;
    data_turnoff.Vds = Vds_turnoff;
    data_turnoff.Id = Id_turnoff;
    
    data_turnon = table();
    data_turnon.t = t_turnon;
    data_turnon.Vgs = Vgs_turnon;
    data_turnon.Vds = Vds_turnon;
    data_turnon.Id = Id_turnon;
    
    data_trs = table();
    data_trs.turnoff = data_turnoff;
    data_trs.turnon = data_turnon;
    
    data_raw = table();
    data_raw.t = t;
    data_raw.Vgs = Vgs;
    data_raw.Vds = Vds;
    data_raw.Id = Id;

    clear t_turnoff Vgs_turnoff Vds_turnoff Id_turnoff t_turnon Vgs_turnon Vds_turnon Id_turnon
    clear t Vgs Vds Id
    clear data_turnoff data_turnon

    % Perform FFT
    t = data_raw.t;
    Vgs_turnoff = data_trs.turnoff.Vgs;
    Vds_turnoff = data_trs.turnoff.Vds;
    Id_turnoff = data_trs.turnoff.Id;
    Vgs_turnon = data_trs.turnon.Vgs;
    Vds_turnon = data_trs.turnon.Vds;
    Id_turnon = data_trs.turnon.Id;

    fs = length(t)/(t(length(t)) - t(1));
    
    % For Turn off Transient
    N_fft = length(Vgs_turnoff);
    fft_Id_turnoff = fft(Id_turnoff, N_fft);
    fft_Vds_turnoff = fft(Vds_turnoff, N_fft);
    fft_Vgs_turnoff = fft(Vgs_turnoff, N_fft);

    % For Turn on Transient
    N_fft = length(Vgs_turnon);
    fft_Id_turnon = fft(Id_turnon, N_fft);
    fft_Vds_turnon = fft(Vds_turnon, N_fft);
    fft_Vgs_turnon = fft(Vgs_turnon, N_fft);
    
    % Calculate the corresponding frequency vector

    N_turnoff = N_fft;
    f_turnoff = ((0:(N_turnoff-1))*(fs/N_turnoff))';
    N_turnon = N_fft;
    f_turnon = ((0:(N_turnon-1))*(fs/N_turnon))';
    clear t
    clear Vgs_turnoff Vds_turnoff Id_turnoff
    clear Vgs_turnon Vds_turnon Id_turnon
    clear N_turnoff N_turnon N_fft fs
    
    fft_turnoff = table();
    fft_turnoff.f = f_turnoff;
    fft_turnoff.Vgs = fft_Vgs_turnoff;
    fft_turnoff.Vds = fft_Vds_turnoff;
    fft_turnoff.Id = fft_Id_turnoff;
    
    fft_turnon = table();
    fft_turnon.f = f_turnon;
    fft_turnon.Vgs = fft_Vgs_turnon;
    fft_turnon.Vds = fft_Vds_turnon;
    fft_turnon.Id = fft_Id_turnon;
    
    fft_trs = table();
    fft_trs.turnoff = fft_turnoff;
    fft_trs.turnon = fft_turnon;

    clear f_turnoff_zm fft_Vgs_turnoff_zm fft_Vds_turnoff_zm fft_Id_turnoff_zm
    clear f_turnon_zm fft_Vgs_turnon_zm fft_Vds_turnon_zm fft_Id_turnon_zm
    clear fft_turnoff_zm fft_turnon_zm
    
    % Turn-off Data
    abs_fft_Id_turnoff = abs(fft_trs_zm.turnoff.Id)';
    abs_fft_Vds_turnoff = abs(fft_trs_zm.turnoff.Vds)';
    abs_fft_Vgs_turnoff = abs(fft_trs_zm.turnoff.Vgs)';
    
    angle_fft_Id_turnoff = angle(fft_trs_zm.turnoff.Id)';
    angle_fft_Vds_turnoff = angle(fft_trs_zm.turnoff.Vds)';
    angle_fft_Vgs_turnoff = angle(fft_trs_zm.turnoff.Vgs)';
    
    Result.turnoff = [abs_fft_Id_turnoff; abs_fft_Vds_turnoff; abs_fft_Vgs_turnoff;...
                     angle_fft_Id_turnoff; angle_fft_Vds_turnoff; angle_fft_Vgs_turnoff];
    
    
    % Turn-on Data
    abs_fft_Id_turnon = abs(fft_trs_zm.turnon.Id)';
    abs_fft_Vds_turnon = abs(fft_trs_zm.turnon.Vds)';
    abs_fft_Vgs_turnon = abs(fft_trs_zm.turnon.Vgs)';
    
    angle_fft_Id_turnon = angle(fft_trs_zm.turnon.Id)';
    angle_fft_Vds_turnon = angle(fft_trs_zm.turnon.Vds)';
    angle_fft_Vgs_turnon = angle(fft_trs_zm.turnon.Vgs)';
    
    Result.turnon = [abs_fft_Id_turnon; abs_fft_Vds_turnon; abs_fft_Vgs_turnon;...
                     angle_fft_Id_turnon; angle_fft_Vds_turnon; angle_fft_Vgs_turnon];

    % SAVE CSV file

    Vgs_on = 0.1*str2double(name_parts{6});
    Vgs_off = 0.1*str2double(name_parts{8});

    Save_name= [Result_Path,'/Turn_off/','Ron_',num2str(R_on),'_Roff_',num2str(R_off),'_Pulse_',num2str(name_parts{2}),...
        '_Vds_',num2str(name_parts{4}),'_Vgson_',num2str(Vgs_on),'_Vgsoff_',num2str(Vgs_off),...
        '_Resamplefac_',num2str(resample_factor),'_id_',num2str(i-2),'offFFT.csv'];
    
    writematrix(Result.turnoff,Save_name)

    Save_name= [Result_Path,'/Turn_on/','Ron_',num2str(R_on),'_Roff_',num2str(R_off),'_Pulse_',num2str(name_parts{2}),...
    '_Vds_',num2str(name_parts{4}),'_Vgson_',num2str(Vgs_on),'_Vgsoff_',num2str(Vgs_off),...
    '_Resamplefac_',num2str(resample_factor),'_id_',num2str(i-2),'onFFT.csv'];

    writematrix(Result.turnon,Save_name)
    
    
    proctime = toc(iterTimer);
    totalTime = toc(totalTimer);
    formatted_time = datestr(totalTime / (24 * 60 * 60), 'HH:MM:SS');
    fprintf('Iteration %d / %d is done. Process time: %.3f Elapse Time: %s\n', i-2, Data_num-2, proctime,formatted_time);
end