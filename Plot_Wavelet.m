function Plot_Wavelet(t, t_limit, y, cfs, f)
    figure;
    subplot(3,1,1);
    plot(t, y);
    title("Original Data");
    xlim(t_limit);
    
    subplot(3,1,2);
    imagesc([t(1) t(length(t))], [f(1) f(length(f)) ], abs(cfs)); % Plot the CWT coefficients using imagesc
    colormap(jet); % Choose the colormap for the CWT coefficients plot
    % colorbar; % Show the color bar for the CWT coefficients
    title("CWT abs");
    xlim(t_limit);

    subplot(3,1,3);
    imagesc([t(1) t(length(t))], [f(1) f(length(f)) ], angle(cfs)); % Plot the CWT coefficients using imagesc
    colormap(jet); % Choose the colormap for the CWT coefficients plot
    % colorbar; % Show the color bar for the CWT coefficients
    title("CWT angle");
    xlim(t_limit);