function  View_Dispersion(w_0_n, gamma_n,sigma_n,eps_infin,lambda_min,lambda_max,num_points,material)
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
    

    eps_mag = (eps_re.^2 + eps_im.^2).^(1/2); % Magnitude of eps
    n_0 = ((eps_mag + eps_re)./2).^(1/2);
    kappa = ((eps_mag - eps_re)./2).^(1/2);
    
    title_1 = material + " n_0 by Wavelength";
    figure(1)
    plot(lambda_vec,n_0)
    title(title_1)
    xlabel('lambda (um)')
    ylabel('n_0')

    title_2 = material + " n_0 by Energy";
    figure(2)
    plot(Evs,n_0)
    title(title_2)
    xlabel('E (eV)')
    ylabel('n_0')

    title_3 = material + " Real epsilon by Wavelength";
    figure(3)
    plot(lambda_vec,eps_re)
    title(title_3)
    xlabel('lambda (um)')
    ylabel('real eps')

    title_4 = material + " Imaginary Epsilon by Wavelength";
    figure(4)
    plot(lambda_vec,eps_im)
    title(title_4)
    xlabel('lambda (um)')
    ylabel('im eps')

    title_5 = material + " Extintion Coefficient by Wavelength";
    figure(5)
    plot(lambda_vec,kappa)
    title(title_5)
    xlabel('lambda (um)')
    ylabel('kappa')

    title_6 = material + " Extinction Coefficient by Energy";
    figure(6)
    plot(Evs,kappa)
    title(title_6)
    xlabel('E (eV')
    ylabel('kappa')

else
    warning("Error: Make sure that the length of w_0_n = gamma_n = sigma_n \n \n")
end

end
