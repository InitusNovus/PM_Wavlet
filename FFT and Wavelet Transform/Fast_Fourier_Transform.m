%% TODO:
%% Import CSV File
clearvars;
% Write HyperParameter
File_PATH = '../Data/Ron_50_Roff_50';
Result_Path = './../Data/FFT_Data';
Select_variable='Test';

% Please write Resist
R_on=50;
R_off=50;

File_List = dir(File_PATH);
[Data_num,~] = size(File_List);

% Plot figure
% fig_handle = figure;
fig_handle = figure(1);
%% File Selection Section
totalTimer = tic; % Total elapsed timer
for k = 3: Data_num % The first and the second file is "." and "..".
% for k = 152   %For test
    %% Process Started
    iterTimer = tic; % Iteration timer
    %% File Parsing

    % File name configuration
    file_name = File_List(k).name;
    name_parts = split(file_name,'_');
    file_name = strcat(File_PATH,'/',file_name);
    data = readtable(file_name);

    % Data parsing
    data = renamevars(data, "Labels", "Time");
    t = data.Time;
    Vgs = data.Vgs;
    Vds = data.Vds;
    Id = data.Id;
    %% Finding the Switching Transients

    % Eliminate the mean value in Vds to find the switching transient
    Vds_mean = mean(Vds);
    Vds_zeromean = Vds - Vds_mean;
    min_gap = 10000; % Arbitrary value

    zero_crossings = Find_zeroCrossing(Vds_zeromean, min_gap);
    assert(length(zero_crossings) == 4, "Error: Switching transient detection failed");

    % Interpolate the time values corresponding to the zero crossings
    zero_crossing_times = interp1(1:length(t), t, zero_crossings);
    clear zero_crossings

    zero_crossing_names = {"Charging"; "Turn-Off"; "Turn-On"; "End"};
    st_time = containers.Map();
    for i = 1:4
        st_time(zero_crossing_names{i}) = zero_crossing_times(i);
    end
    clear i
    
    % clear
    clear zero_crossing_times zero_crossing_names
    clear Vds_mean Vds_zeromean
    clear min_gap

    % Find the switching transients
    fs = length(t)/(t(length(t)) - t(1));
    testTime = 1e-6; % seconds
    testLength = floor(fs * testTime + eps);
    % st_time_range = 0.1 * [-1e-6 9e-6];
    % st_time_range = 0.1 * [(-0.3e-6) (-0.3e-6 + 10e-6)];
    % st_time_range = 0.1 * [1e-6 11e-6];
    % st_time_range = 0.1 * [0.5e-6 10.5e-6];

    st_time_range_turnoff = 0.1 * [-1e-6 9e-6];
    % st_time_range_turnon = 0.1 * [1e-6 11e-6];
    st_time_range_turnon = 0.1 * [-1e-6 9e-6];

    st_time_range_turnoff = st_time_range_turnoff + eps;
    st_time_range_turnon = st_time_range_turnon + eps;

    % Turn-Off transient
    turnoff_time_search_range = st_time("Turn-Off") + st_time_range_turnoff;
    index_search_turnoff = (turnoff_time_search_range(1) < t) & (t < turnoff_time_search_range(2));

    Id_search_turnoff = Id(index_search_turnoff);
    min_gap = 10000; % Arbitrary value
    t_search_turnoff_zcIndex = floor(Find_zeroCrossing(Id_search_turnoff, min_gap));
    t_search_turnoff_startIndex = find(index_search_turnoff, 1);
    t_turnoff_startIndex = t_search_turnoff_startIndex + t_search_turnoff_zcIndex;
    index_turnoff = t_turnoff_startIndex + (0:(testLength-1));

    t_turnoff = t(index_turnoff);
    Vgs_turnoff = Vgs(index_turnoff);
    Vds_turnoff = Vds(index_turnoff);
    Id_turnoff = Id(index_turnoff);
    
    % clear
    clear turnoff_time_search_range index_search_turnoff
    clear t_search_turnoff Id_search_turnoff min_gap
    clear t_search_turnoff_zcIndex t_search_turnoff_startIndex t_turnoff_startIndex
    clear index_turnoff

    % Turn-On transient
    turnon_time_search_range = st_time("Turn-On") + st_time_range_turnon;
    index_search_turnon = (turnon_time_search_range(1) < t) & (t < turnon_time_search_range(2));

    Vds_search_turnon = Vds(index_search_turnon);
    min_gap = 10000; % Arbitrary value
    t_search_turnon_zcIndex = floor(Find_zeroCrossing(Vds_search_turnon, min_gap));

    t_search_turnon_startIndex = find(index_search_turnon, 1);
    t_turnon_startIndex = t_search_turnon_startIndex + t_search_turnon_zcIndex;
    index_turnon = t_turnon_startIndex + (0:(testLength-1));

    t_turnon = t(index_turnon);
    Vgs_turnon = Vgs(index_turnon);
    Vds_turnon = Vds(index_turnon);
    Id_turnon = Id(index_turnon);
    
    % clear
    clear turnon_time_search_range index_search_turnon
    clear t_search_turnon Vds_search_turnon min_gap
    clear t_search_turnon_zcIndex t_search_turnon_startIndex t_turnon_startIndex
    clear index_turnon

    % clear

    clear turnoff_time_range turnon_time_range
    clear index_turnoff index_turnon
    clear st_time_range

    %% Data to Table

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
    
    % clear
    clear t_turnoff Vgs_turnoff Vds_turnoff Id_turnoff t_turnon Vgs_turnon Vds_turnon Id_turnon
    clear t Vgs Vds Id
    clear data_turnoff data_turnon

    %% Zeromean
    data_turnoff_arr_with_t = table2array(data_trs.turnoff);
    data_turnoff_arr = data_turnoff_arr_with_t(:, 2:end);
    clear data_turnoff_arr_with_t
    data_turnoff_mean = mean(data_turnoff_arr);
    data_turnoff_zm = data_turnoff_arr - data_turnoff_mean;

    data_trs_turnoff_zm = data_trs.turnoff;
    data_trs_turnoff_zm{:, 2:4} = data_turnoff_zm;

    data_trs_zm = table();
    data_trs_zm.turnoff = data_trs_turnoff_zm;

    data_turnon_arr_with_t = table2array(data_trs.turnon);
    data_turnon_arr = data_turnon_arr_with_t(:, 2:end);
    clear data_turnon_arr_with_t
    data_turnon_mean = mean(data_turnon_arr);
    data_turnon_zm = data_turnon_arr - data_turnon_mean;

    data_trs_turnon_zm = data_trs.turnon;
    data_trs_turnon_zm{:, 2:4} = data_turnon_zm;

    data_trs_zm.turnon = data_trs_turnon_zm;
    
    %clear
    clear data_turnoff_arr data_turnoff_mean data_turnoff_zm data_trs_turnoff_zm
    clear data_turnon_arr data_turnon_mean data_turnon_zm data_trs_turnon_zm

    %% FFT (Odler version)

    % % Parse Internal Variables
    % t = data_raw.t;
    % Vgs_turnoff = data_trs_zm.turnoff.Vgs;
    % Vds_turnoff = data_trs_zm.turnoff.Vds;
    % Id_turnoff = data_trs_zm.turnoff.Id;
    % Vgs_turnon = data_trs_zm.turnon.Vgs;
    % Vds_turnon = data_trs_zm.turnon.Vds;
    % Id_turnon = data_trs_zm.turnon.Id;
    % 
    % % Sampling Frequency
    % fs = length(t)/(t(length(t)) - t(1));
    % 
    % % For Turn off Transient
    % N_fft = length(Vgs_turnoff);
    % fft_Id_turnoff_zm = fft(Id_turnoff, N_fft);
    % fft_Vds_turnoff_zm = fft(Vds_turnoff, N_fft);
    % fft_Vgs_turnoff_zm = fft(Vgs_turnoff, N_fft);
    % 
    % % For Turn on Transient
    % 
    % N_fft = length(Vgs_turnon);
    % fft_Id_turnon_zm = fft(Id_turnon, N_fft);
    % fft_Vds_turnon_zm = fft(Vds_turnon, N_fft);
    % fft_Vgs_turnon_zm = fft(Vgs_turnon, N_fft);
    % 
    % % Calculate the corresponding grequency vector
    % 
    % N_turnoff = N_fft;
    % f_turnoff_zm = ((0:(N_turnoff-1))*(fs/N_turnoff))';
    % N_turnon = N_fft;
    % f_turnon_zm = ((0:(N_turnon-1))*(fs/N_turnon))';
    % clear t
    % clear Vgs_turnoff Vds_turnoff Id_turnoff
    % clear Vgs_turnon Vds_turnon Id_turnon
    % clear N_turnoff N_turnon N_fft fs
    % 
    % fft_turnoff_zm = table();
    % fft_turnoff_zm.f = f_turnoff_zm;
    % fft_turnoff_zm.Vgs = fft_Vgs_turnoff_zm;
    % fft_turnoff_zm.Vds = fft_Vds_turnoff_zm;
    % fft_turnoff_zm.Id = fft_Id_turnoff_zm;
    % 
    % fft_turnon_zm = table();
    % fft_turnon_zm.f = f_turnon_zm;
    % fft_turnon_zm.Vgs = fft_Vgs_turnon_zm;
    % fft_turnon_zm.Vds = fft_Vds_turnon_zm;
    % fft_turnon_zm.Id = fft_Id_turnon_zm;
    % 
    % fft_trs_zm = table();
    % fft_trs_zm.turnoff = fft_turnoff_zm;
    % fft_trs_zm.turnon = fft_turnon_zm;
    % clear f_turnoff_zm fft_Vgs_turnoff_zm fft_Vds_turnoff_zm fft_Id_turnoff_zm
    % clear f_turnon_zm fft_Vgs_turnon_zm fft_Vds_turnon_zm fft_Id_turnon_zm
    % clear fft_turnoff_zm fft_turnon_zm
    % 
    % % Turn-on Data
    % abs_fft_Id_turnon = abs(fft_trs_zm.turnon.Id)';
    % abs_fft_Vds_turnon = abs(fft_trs_zm.turnon.Vds)';
    % abs_fft_Vgs_turnon = abs(fft_trs_zm.turnon.Vgs)';
    % 
    % angle_fft_Id_turnon = angle(fft_trs_zm.turnon.Id)';
    % angle_fft_Vds_turnon = angle(fft_trs_zm.turnon.Vds)';
    % angle_fft_Vgs_turnon = angle(fft_trs_zm.turnon.Vgs)';
    % 
    % Result.turnon = [abs_fft_Id_turnon; abs_fft_Vds_turnon; abs_fft_Vgs_turnon;...
    %     angle_fft_Id_turnon; angle_fft_Vds_turnon; angle_fft_Vgs_turnon];
    % 
    % 
    % % Turn-off Data
    % abs_fft_Id_turnoff = abs(fft_trs_zm.turnoff.Id)';
    % abs_fft_Vds_turnoff = abs(fft_trs_zm.turnoff.Vds)';
    % abs_fft_Vgs_turnoff = abs(fft_trs_zm.turnoff.Vgs)';
    % 
    % angle_fft_Id_turnoff = angle(fft_trs_zm.turnoff.Id)';
    % angle_fft_Vds_turnoff = angle(fft_trs_zm.turnoff.Vds)';
    % angle_fft_Vgs_turnoff = angle(fft_trs_zm.turnoff.Vgs)';
    % 
    % Result.turnoff = [abs_fft_Id_turnoff; abs_fft_Vds_turnoff; abs_fft_Vgs_turnoff;...
    %     angle_fft_Id_turnoff; angle_fft_Vds_turnoff; angle_fft_Vgs_turnoff];


    %% Welch Method

    % Basic Parameters
    fs = length(data.Time)/(data.Time(end) - data.Time(1));
    trs = ["turnoff", "turnon"];
    ch = ["Id", "Vds", "Vgs"];

    % [Hyperparameters] for pwelch
    % "winTime = dataTime" and same Data length: 
    %       Just FFT with Hanning window and normalized amplitude squared
    winTime = 1000; %ns
    overlabRatio = 0.25; %ratio
    nfftRatio = 10; %ratio

    for whatTrs = trs
        data_internal = data_trs_zm.(whatTrs);
        dataTime = (length(data_internal.t)/fs) * 1e9; %ns
        % disp("dataTime = " + dataTime + " ns"); disp("winTime = " + winTime + " ns");
        win = round(winTime*1e-9 * fs);
        nfft = ceil(nfftRatio * win);
        noverlap = ceil(overlabRatio * win);
        assert(win > 0);

        for ch_name = ch
            [psd.(whatTrs).(ch_name).pxx, psd.(whatTrs).(ch_name).f, psd.(whatTrs).(ch_name).pxxc] = ...
                pwelch(single(data_internal.(ch_name)), win, noverlap, nfft, fs, 'onesided', 'ConfidenceLevel', 0.95);
        end
    end
    %% Process Finished
    proctime = toc(iterTimer);
    %% Writing Started
    writeTimer = tic;
    %% Result Parsing
    % (Id; Vds; Vgs) Each data for channel is in the rows.
    clear write
    write_range = 21:1001;
    write.turnoff(1,:) = psd.turnoff.Id.pxx(write_range)';
    write.turnoff(2,:) = psd.turnoff.Vds.pxx(write_range)';
    write.turnoff(3,:) = psd.turnoff.Vgs.pxx(write_range)';
    write.turnon(1,:) = psd.turnon.Id.pxx(write_range)';
    write.turnon(2,:) = psd.turnon.Vds.pxx(write_range)';
    write.turnon(3,:) = psd.turnon.Vgs.pxx(write_range)';
    %% Plot for the Test
    fig_handle;
    clf;
    for j = 1:length(trs)
        whatTrs = trs(j);
        % figure(i);
        for i = 1:3
            ch_name = ch(i);

            subplot(3, length(trs), 2*i - j + 1);
            hold on
            plot(psd.(whatTrs).(ch_name).f(write_range) * 1e-6, write.(whatTrs)(i,:), LineWidth=1);
            hold off
            grid on;
            axis('tight')
            ylabel(ch_name)
        end
        subplot(3, length(trs), 2*1 - j + 1);
        title(whatTrs);
        subplot(3, length(trs), 2*3 - j + 1);
        xlabel("Frequency [MHz]")

        drawnow;
    end
    %% SAVE CSV file

    Vgs_on = 0.1*str2double(name_parts{6});
    Vgs_off = 0.1*str2double(name_parts{8});
    data_id = k - 2;
    


    Save_name= [Result_Path,'/Turn_off/','Ron_',num2str(R_on),'_Roff_',num2str(R_off),'_Pulse_',...
        num2str(name_parts{2}),'_Vds_',num2str(name_parts{4}),'_Vgson_',num2str(Vgs_on),'_Vgsoff_',num2str(Vgs_off),...
        '_Turnoff_',num2str(0),'_id_',num2str(data_id),'.FFT.csv'];

    writematrix(write.turnoff,Save_name)

    Save_name= [Result_Path,'/Turn_on/','Ron_',num2str(R_on),'_Roff_',num2str(R_off),'_Pulse_',...
        num2str(name_parts{2}),'_Vds_',num2str(name_parts{4}),'_Vgson_',num2str(Vgs_on),'_Vgsoff_',num2str(Vgs_off),...
        '_Turnon_',num2str(1),'_id_',num2str(data_id),'.FFT.csv'];

    writematrix(write.turnon,Save_name)


    writeTime = toc(writeTimer);
    totalTime = toc(totalTimer);
    formatted_time = datestr(totalTime / (24 * 60 * 60), 'HH:MM:SS');
    fprintf('Iteration %d / %d is done.\n-Process time: %.3f Write time: %.3f\n-Total Elapse Time: %s\n',...
        data_id, Data_num-2, proctime, writeTime, formatted_time);
end