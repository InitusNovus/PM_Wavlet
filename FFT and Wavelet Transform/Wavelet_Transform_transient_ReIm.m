%% TODO:
%% Import CSV File
clearvars;
disp("[Wavelet Transform Data Processor]")
% Write HyperParameter
File_PATH = "../Data/Ron_50_Roff_50";
Result_Path = "./../Data/Wavelet_Result_Transient";
Select_variable="Test";
resFactor = 50;

writingOriginalTrs = true;
writingTrsAbsAngle = true;
writingTrsRealImag = true;

folderOriginalTrs = "/TimeDomain";
folderTrsAbsAngle = "/Wavelet";
folderTrsRealImag = "/Wavelet_ReIm";

% Please write Resist
R_on=50;
R_off=50;

File_List = dir(File_PATH);
[Data_num,~] = size(File_List);

% Plot figure
% fig_handle = figure;
% fig_handle = figure(1);
%% 
fprintf("- Source Folder: %s\n", File_PATH);
fprintf("- Destination Folder: %s\n", Result_Path);
disp("Outputs: ")
if writingOriginalTrs
disp("- Time Domain Transient Waveform")
fprintf("   at %s\n", folderOriginalTrs)
end
if writingTrsAbsAngle
disp("- Wavelet Magnitude and Phase angle")
fprintf("   at %s\n", folderTrsAbsAngle)
end
if writingTrsRealImag
disp("- Wavelet Real and Imag parts")
fprintf("   at %s\n", folderTrsRealImag)
end
%% File Selection Section
disp("Process Starting");
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
    %% Find the switching transient
    % Eliminate the mean value in Vds to find the switching transient by means of finding the zerocrossing of Vds waveform.
    Vds_mean = mean(Vds);
    Vds_zeromean = Vds - Vds_mean;

    % Find the zerocrossing
    min_gap = 10000; % Arbitrary value

    zero_crossings = Find_zeroCrossing(Vds_zeromean, min_gap);
    assert(length(zero_crossings) == 4, "Error: Switching transient detection failed")

    % Interpolate the time values corresponding to the zero crossings
    zero_crossing_times = interp1(1:length(t), t, zero_crossings);
    clear zero_crossings
    
    % Parse the zerocrossing point to switching transient
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
    
    %% Find the switching transients
    percentOffset_turnoff = -25.0;
    percentOffset_turnon = -25.0;

    % testTime = 1.392e-6; % seconds  % To be test length 8700
    testTime = 2e-6; % seconds  % To be test length 10000
    % testTime = 4.0e-6; % seconds  %

    % st_time_range_turnoff_offset = -0.35e-6;
    % st_time_range_turnon_offset = -0.35e-6;
    st_time_range_turnoff_offset = testTime * percentOffset_turnoff/100;
    st_time_range_turnon_offset = testTime * percentOffset_turnon/100;
    fs = length(t)/(t(length(t)) - t(1));
    testLength = floor(fs * testTime + eps);

    st_time_range_turnoff = [0 testTime] + st_time_range_turnoff_offset;
    st_time_range_turnon = [0 testTime] + st_time_range_turnon_offset;

    st_time_range_turnoff = st_time_range_turnoff + eps;
    st_time_range_turnon = st_time_range_turnon + eps;
    % Turn-Off transient
    turnoff_time_search_range = st_time("Turn-Off") + st_time_range_turnoff;
    index_search_turnoff = (turnoff_time_search_range(1) < t) & (t < turnoff_time_search_range(2));
    index_turnoff = index_search_turnoff;

    t_turnoff = t(index_turnoff);
    Vgs_turnoff = Vgs(index_turnoff);
    Vds_turnoff = Vds(index_turnoff);
    Id_turnoff = Id(index_turnoff);
    
    %clear: Turn-Off transient
    clear turnoff_time_search_range index_search_turnoff
    clear t_search_turnoff Id_search_turnoff min_gap
    clear t_search_turnoff_zcIndex t_search_turnoff_startIndex t_turnoff_startIndex
    clear percentOffset_turnoff st_time_range_turnoff st_time_range_turnoff_offset

    % Turn-On transient
    turnon_time_search_range = st_time("Turn-On") + st_time_range_turnon;
    index_search_turnon = (turnon_time_search_range(1) < t) & (t < turnon_time_search_range(2));
    index_turnon = index_search_turnon;

    t_turnon = t(index_turnon);
    Vgs_turnon = Vgs(index_turnon);
    Vds_turnon = Vds(index_turnon);
    Id_turnon = Id(index_turnon);

    % clear: Turn-On transient
    clear turnon_time_search_range index_search_turnon
    clear t_search_turnon Vds_search_turnon min_gap
    clear t_search_turnon_zcIndex t_search_turnon_startIndex t_turnon_startIndex
    clear percentOffset_turnon st_time_range_turnon st_time_range_turnon_offset

    % clear
    clear turnoff_time_search_range turnon_time_range
    clear st_time_range testTime testLength st_time
    
    %% Data to table
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

    data_trs.turnoff = data_turnoff;
    data_trs.turnon = data_turnon;
    data_trs.fs = fs;
    data_trs.index.turnoff = index_turnoff;
    data_trs.index.turnon = index_turnon;

    data_raw.data.t = t;
    data_raw.data.Vgs = Vgs;
    data_raw.data.Vds = Vds;
    data_raw.data.Id = Id;
    data_raw.fs = fs;
    
    % clear
    clear t_turnoff Vgs_turnoff Vds_turnoff Id_turnoff t_turnon Vgs_turnon Vds_turnon Id_turnon
    clear t Vgs Vds Id
    clear data_turnoff data_turnon
    clear index_turnoff index_turnon

    %% Continuous Wavelet Transform

    % Parameter Setting
    trs = ["turnoff" "turnon"];
    ch = ["Id" "Vds" "Vgs"];
    data_trs.cwt = struct;
    data_trs.cwt_res = struct;

    for i = 1:length(trs)
        whatTrs = trs(i);
        % Load local variables
        t = data_trs.(whatTrs).t;
        fs = data_trs.fs;
        idx = data_trs.index.(whatTrs);

        % For each channel
        for j = 1:length(ch)
            whatCh = ch(j);
            waveform = data_trs.(whatTrs).(whatCh);
            % For entire waveform
            [cfs_full, f] = cwt(data_raw.data.(whatCh), 'Morse', fs, VoicePerOctave = 5);

            % Cut for the transient range
            cfs = single(cfs_full(:, idx));
            data_trs.cwt.(whatTrs).(whatCh).cfs = cfs;
            data_trs.cwt.(whatTrs).(whatCh).f = f;

            % Resampling
            cfs_res = resample(cfs.', 1, resFactor).';
            data_trs.cwt_res.resFactor = resFactor;
            data_trs.cwt_res.(whatTrs).(whatCh).cfs = cfs_res;
            data_trs.cwt_res.(whatTrs).(whatCh).f = f;
        end
    end
    clear i j
    clear t fs idx
    clear whatTra whatCh cfs_full cfs cfs_res

    %% Process Finished
    proctime = toc(iterTimer);
    %% Writing Started
    writeTimer = tic;
    %% Result Parsing
    clear write
    clear write_timeTrs
    clear write_ReIm
    
    % Time Domain Transient
    for whatTrs = trs
        write_timeTrs.(whatTrs) = table();
        for whatCh = ["t" ch]
            write_timeTrs.(whatTrs).(whatCh) = data_trs.(whatTrs).(whatCh);
        end
    end

    % Wavelet Transient: 
    % vertcat(abs(Id), abs(Vds), abs(Vgs), angle(Id), angle(Vds), angle(Vgs))
    for whatTrs = trs
        write.(whatTrs) = [];
        for whatCh = ch
            write.(whatTrs) = [write.(whatTrs); abs(data_trs.cwt_res.(whatTrs).(whatCh).cfs)];
        end
        for whatCh = ch
            write.(whatTrs) = [write.(whatTrs); angle(data_trs.cwt_res.(whatTrs).(whatCh).cfs)];
        end
    end

    % Wavelet Transient Real and Imag: 
    % vertcat(real(Id), real(Vds), real(Vgs), imag(Id), imag(Vds), imag(Vgs))
    for whatTrs = trs
        write_ReIm.(whatTrs) = [];
        for whatCh = ch
            write_ReIm.(whatTrs) = [write_ReIm.(whatTrs); real(data_trs.cwt_res.(whatTrs).(whatCh).cfs)];
        end
        for whatCh = ch
            write_ReIm.(whatTrs) = [write_ReIm.(whatTrs); imag(data_trs.cwt_res.(whatTrs).(whatCh).cfs)];
        end
    end

    % %% Plot for the Test
    % fig_handle;
    % clf;
    % for j = 1:length(trs)
    %     whatTrs = trs(j);
    %     % figure(i);
    %     for i = 1:3
    %         ch_name = ch(i);
    % 
    %         subplot(3, length(trs), 2*i - j + 1);
    %         hold on
    %         plot(write_timeTrs.(whatTrs).t, write_timeTrs.(whatTrs).(ch_name, LineWidth=1);
    %         hold off
    %         grid on;
    %         axis('tight')
    %         ylabel(ch_name)
    %     end
    %     subplot(3, length(trs), 2*1 - j + 1);
    %     title(whatTrs);
    %     subplot(3, length(trs), 2*3 - j + 1);
    %     xlabel("Frequency [MHz]")
    % 
    %     drawnow;
    % end
    %%  Save CSV file
    Vgs_on = 0.1*str2double(name_parts{6});
    Vgs_off = 0.1*str2double(name_parts{8});
    data_id = k - 2;
    
    % Time Domain Transient
    if writingOriginalTrs
        for whatTrs = trs
            Save_name= Result_Path + folderOriginalTrs + "/" + whatTrs + "/" + "Ron_" + num2str(R_on) + "_Roff_" + num2str(R_off) + "_Pulse_" + ...
                num2str(name_parts{2}) + "_Vds_" + num2str(name_parts{4}) + "_Vgson_" + num2str(Vgs_on) + ...
                "_Vgsoff_" + num2str(Vgs_off) +  "_" + whatTrs + "_" + num2str(0) + "_id_" + num2str(data_id) + ".TRS.csv";
            writetable(write_timeTrs.(whatTrs), Save_name)
        end
    end

    % Wavelet Transient
    if writingTrsAbsAngle
        for whatTrs = trs
            Save_name= Result_Path + folderTrsAbsAngle + "/" + whatTrs + "/" + "Ron_" + num2str(R_on) + "_Roff_" + num2str(R_off) + "_Pulse_" + ...
                num2str(name_parts{2}) + "_Vds_" + num2str(name_parts{4}) + "_Vgson_" + num2str(Vgs_on) + ...
                "_Vgsoff_" + num2str(Vgs_off) +  "_" + whatTrs + "_" + num2str(0) + "_id_" + num2str(data_id) + ".WVL.csv";
            writematrix(write.(whatTrs), Save_name)
        end
    end
    
    % Real parts and Imaginary parts
    if writingTrsRealImag
        for whatTrs = trs
            Save_name= Result_Path + folderTrsRealImag + "/" + whatTrs + "/" + "Ron_" + num2str(R_on) + "_Roff_" + num2str(R_off) + "_Pulse_" + ...
                num2str(name_parts{2}) + "_Vds_" + num2str(name_parts{4}) + "_Vgson_" + num2str(Vgs_on) + ...
                "_Vgsoff_" + num2str(Vgs_off) +  "_" + whatTrs + "_" + num2str(0) + "_id_" + num2str(data_id) + ".WVL_RI.csv";
            writematrix(write_ReIm.(whatTrs), Save_name)
        end
    end
    writeTime = toc(writeTimer);
    %% Disp Elapsed Time
    totalTime = toc(totalTimer);
    formatted_time = datestr(totalTime / (24 * 60 * 60), 'HH:MM:SS');
    fprintf('Iteration %d / %d is done.\n- Process time: %.3f Write time: %.3f\n- Total Elapse Time: %s\n',...
        data_id, Data_num-2, proctime, writeTime, formatted_time);
end
























