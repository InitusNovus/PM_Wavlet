%% Import CSV

% Write HyperParameter
File_PATH_GAN = '../Data/Gan_Data/FFT/Result.csv';
data_Gan = readtable(File_PATH_GAN);

File_PATH_Original = '../Data/Gan_Data/FFT/Original.csv';
data_Original = readtable(File_PATH_Original);


batch = size(data_Gan,1)/3;
data_Gan = reshape(data_Gan{:,:},batch,[],3);
data_Original = reshape(data_Original{:,:},batch,3,[]);

check_data=data_Gan(1,:,:);
check_data =reshape(check_data,981,3);

for i = 1:3
    subplot(2,3,i)
    plot(reshape(data_Gan(1,i,:),1,[]))
    title('Generate_data')

    subplot(2,3,3+i)

    plot(reshape(data_Original(1,i,:),1,[]))
    title('Original')

end
