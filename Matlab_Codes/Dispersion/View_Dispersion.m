function  View_Dispersion(w_0_n, gamma_n,sigma_n,eps_infin,lambda_min,lambda_max,num_points,material,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function plots a dispersion relation as it will be calculated by,
%   and implemented in, meep.
%   
%   Length scale is set to 1 um
%
%   w_0_n -> A vector of the resonant frequncies of the model (in meep 
%   units 1 um / lambda um)
%   
%   gamma_n -> A vector of the gamma parameter needed for a lorentz model
%
%   sigma_n -> A vector of sigma values used in scaling each lorentz
%   resonance
%
%   eps_infin -> Epsilon value at inifinty, a common input to a lorentian
%   model
%
%   lambda_min -> Minimum wavelength which model will be calculated for
%
%   lambda_max -> Maximum wavelength which model will be calculated for
%   
%   num_points -> Number of points in the plot
%
%   material -> String giving the name of material you are modeling ex
%   material = "BTO" for BTO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458; % Speed of light m/s
h = 4.1357 * 10^-15; % Planks const eV-s

lambda_vec = transpose(lambda_min:((lambda_max-lambda_min)/num_points):lambda_max); % Wavelength vector
data = data(((1./data(:,1)) >= lambda_min) & ((1./data(:,1)) <= lambda_max),:); % Cutting data to fit lambda range
data_lambda = 1./data(:,1); % Data wavelengths
data_freqs = 1 ./ data_lambda;
data_Evs = h .* c ./ (data_lambda .* 10^-6); % Data energies
freqs = 1 ./ lambda_vec; % Frequency vector (meep units)

Evs = h .* c ./ (lambda_vec .* 10^-6); % Energies vector (eV)

if (length(w_0_n) == length(gamma_n)) && (length(w_0_n) == length(sigma_n))
    eps_re = ones(length(lambda_vec),1) .* eps_infin; % Initial real epsilon vector (set to eps_infin)
    eps_im = zeros(length(lambda_vec),1); % Initial imaginary epsilon vector (initialized to zero)
    for index = 1:1:length(w_0_n)
        eps_re = eps_re + (w_0_n(index).^2 .* sigma_n(index)) .* ...
            (w_0_n(index).^2 - freqs.^2) ./ ((w_0_n(index).^2 - freqs.^2).^2 ...
            + freqs.^2 .* (gamma_n(index) ./ (2 * pi)).^2 );
        eps_im = eps_im + (w_0_n(index).^2 .* sigma_n(index)) .* freqs .*...
            (gamma_n(index) ./ (2 * pi)) ./ ((w_0_n(index).^2 - freqs.^2).^2 ...
            + freqs.^2 .* (gamma_n(index) ./ (2 * pi)).^2 );
    end
    
    data_eps_mag = (data(:,2).^2 + data(:,3).^2).^(1/2); % Magnitude of eps
    data_n_0 = ((data_eps_mag + data(:,2))./2).^(1/2);
    data_kappa = ((data_eps_mag - data(:,2))./2).^(1/2);
    
    eps_mag = (eps_re.^2 + eps_im.^2).^(1/2); % Magnitude of eps
    n_0 = ((eps_mag + eps_re)./2).^(1/2);
    kappa = ((eps_mag - eps_re)./2).^(1/2);
    
    width = 1.5;
    lineColor = [0.15, 0.15, 0.15];
    dataColor = [0.45 0.45 0.45];
    sig_figs = 2;
    data_spacing = 5;
    dataSize = 5;

    idx = 1:data_spacing:length(data_lambda);

    title_1 = material + " Refractive Index by Wavelength";
    figure(1)
    hold on
    plot(lambda_vec,n_0,'Color',lineColor, 'LineWidth',width)
    %scatter(data_lambda(idx), data_n_0(idx), 's', 'Color', dataColor)
    plot(data_lambda(idx), data_n_0(idx), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_1, 'Interpreter','latex')
    xlabel('Wavelength ($\mu$m)', 'Interpreter','latex')
    ylabel('$n_0$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin1, ymax1] = scale_graph(data_n_0, n_0, sig_figs);
    ylim([ymin1, ymax1])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off

    title_2 = material + " Refractive Index by Energy";
    figure(2)
    hold on
    plot(Evs,n_0,'Color',lineColor, 'LineWidth',width)
    %scatter(data_Evs(idx), data_n_0(idx), 's', 'Color', dataColor)
    plot(data_Evs(idx), data_n_0(idx), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_2, 'Interpreter','latex')
    xlabel('Energy (eV)', 'Interpreter','latex')
    ylabel('$n_0$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin2, ymax2] = scale_graph(data_n_0, n_0, sig_figs);
    ylim([ymin2, ymax2])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off

    title_3 = material + " Real $\epsilon$ by Wavelength";
    figure(3)
    hold on
    plot(lambda_vec,eps_re,'Color',lineColor, 'LineWidth',width)
    %scatter(data_lambda(idx), data(idx,2), 's', 'Color',dataColor)
    plot(data_lambda(idx), data(idx,2), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_3, 'Interpreter','latex')
    xlabel('Wavelegnth ($\mu$m)', 'Interpreter','latex')
    ylabel('Real $\epsilon$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin3, ymax3] = scale_graph(data(:,2), eps_re, sig_figs);
    ylim([ymin3, ymax3])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off

    title_4 = material + " Imaginary $\epsilon$ by Wavelength";
    figure(4)
    hold on
    plot(lambda_vec,eps_im,'Color',lineColor, 'LineWidth',width)
    %scatter(data_lambda(idx), data(idx,3), 's', 'Color', dataColor)
    plot(data_lambda(idx), data(idx,3), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_4, 'Interpreter','latex')
    xlabel('Wavelength ($\mu$m)', 'Interpreter','latex')
    ylabel('Imaginary $\epsilon$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin4, ymax4] = scale_graph(data(:,3), eps_im, sig_figs);
    ylim([ymin4, ymax4])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off

    title_5 = material + " Extintion Coefficient by Wavelength";
    figure(5)
    hold on
    plot(lambda_vec,kappa,'Color',lineColor, 'LineWidth',width)
    %scatter(data_lambda(idx), data_kappa(idx), 's', 'Color', dataColor)
    plot(data_lambda(idx), data_kappa(idx), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_5, 'Interpreter','latex')
    xlabel('Wavelength ($\mu$m)', 'Interpreter','latex')
    ylabel('$\kappa$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin5, ymax5] = scale_graph(data_kappa, kappa, sig_figs);
    ylim([ymin5, ymax5])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off

    title_6 = material + " Extinction Coefficient by Energy";
    figure(6)
    hold on
    plot(Evs,kappa,'Color',lineColor, 'LineWidth',width)
    %scatter(data_Evs(idx), data_kappa(idx), 's', 'Color', dataColor)
    plot(data_Evs(idx), data_kappa(idx), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_6, 'Interpreter','latex')
    xlabel('Energy (eV)', 'Interpreter','latex')
    ylabel('$\kappa$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin6, ymax6] = scale_graph(data_kappa, kappa, sig_figs);
    ylim([ymin6, ymax6])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off

    title_7 = material + " Real $\epsilon$ by $\omega$";
    figure(7)
    hold on
    plot(freqs,eps_re,'Color',lineColor, 'LineWidth',width)
    %scatter(data_freqs(idx), data(idx,2), 's', 'Color', dataColor)
    plot(data_freqs(idx), data(idx,2), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_7, 'Interpreter','latex')
    xlabel('$\omega$ (meep units)', 'Interpreter','latex')
    ylabel('Real $\epsilon$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin7, ymax7] = scale_graph(data(:,2), eps_re, sig_figs);
    ylim([ymin7, ymax7])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off

    title_8 = material + " Imaginary $\epsilon$ by $\omega$";
    figure(8)
    hold on
    plot(freqs,eps_im,'Color',lineColor, 'LineWidth',width)
    %scatter(data_freqs, data(:,3), 's', 'Color', dataColor)
    plot(data_freqs(idx), data(idx,3), ...
        's', ...
        'LineStyle', 'none', ...
        'MarkerSize', dataSize, ...
        'MarkerEdgeColor', dataColor, ...
        'MarkerFaceColor', 'none')
    title(title_8, 'Interpreter','latex')
    xlabel('$\omega$ (meep units)', 'Interpreter','latex')
    ylabel('Imaginary $\epsilon$', 'Interpreter','latex')
    legend('Model','Experiment', 'Location', 'best', 'Interpreter','latex')
    legend('boxoff')
    set(gca,'TickLabelInterpreter','latex')
    [ymin8, ymax8] = scale_graph(data(:,3), eps_im, sig_figs);
    ylim([ymin8, ymax8])
    ax = gca;
    ax.LineWidth = 1.2;
    grid on
    hold off


else
    warning("Error: Make sure that the length of w_0_n = gamma_n = sigma_n \n \n")
end

end


function [ymin, ymax] = scale_graph(data, model, sig_figs)
    ymax = max([max(data), max(model)]);
    exponent = floor(log10(abs(ymax)));
    ymax = ymax/10^exponent;
    ymax = ymax*10^(sig_figs - 1);
    ymax = ceil(ymax);
    ymax = ymax/10^(sig_figs - 1);
    ymax = ymax*10^exponent;

    
    ymin = min([min(data), min(model)]);
    exponent = floor(log10(abs(ymin)));
    ymin = ymin/10^exponent;
    ymin = ymin*10^(sig_figs - 1);
    ymin = floor(ymin);
    ymin = ymin/10^(sig_figs - 1);
    ymin = ymin*10^exponent;
end
