function param_opt = fitNlorentzean(data,manual_peak_set,num_peaks, re_frac, im_frac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function currently takes the expiremental data for both real and
% imaginary epsilon plots and the number of peaks desired to fit a
% lorentizan model to the data (using a linear combination of lorenzians).
% it will output the paramters needed to input into meep to recreate this
% model for any given material.
%
% Input:
%
% data -> Experimental data in the form [omega_data, re_eps_data,
% im_eps_data] where omega_data, re_eps_data, and im_eps_data are all
% vectors.
%
% manual_peak_set -> For now THIS MUST BE SET TO 1. This is a boolean value
% to determine if the model will determin the optimal number of peaks
% itself or if it will use the number of peaks given to it be the user. The
% program cannot currently determine this by itself so this must be set to
% 1.
%
% num_peaks -> User determined number of lorenzian peaks wanted in the
% model.
%
% re_frac -> Number (preferably from 1 to 0) that determines the priority
% given to fitting to the real part of the eps data (defualt should be 1)
%
% im_frac -> Number (preferably from 1 to 0) that determines the priority
% given to fitting to the im part of the eps data (default should be 1)
%
% Output:
% 
% param_opt -> The output is of the form: 
%                  [eps_infin;
%                   sigma_1;
%                   omega_1;
%                   gamma_1;
%                   sigma_2;
%                   omega_2;
%                   gamma_2;
%                      :
%                      :
%                   sigma_N;
%                   omega_N;
%                   gamma_N]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_meep_dat = data(:,1);
eps_re_dat = data(:,2);
% re_frac = 1; % How much priority is given to fitting to the real eps data
% im_frac = 1; % How much priority is given to fitting to the im eps data
if(size(data,2)>2)
    eps_im_dat = data(:,3);
end

% Running the following code in the case that the user want a specified
% number of peaks
if(manual_peak_set)
    syms eps_infin real
    syms gamma_n [1 num_peaks]
    assume(gamma_n >= 0);
    syms sigma [1 num_peaks]
    assume(sigma >= 0);
    syms w_0_ [1 num_peaks]
    assume(w_0_ >= 0);
    syms omega
    assume(omega >= 0);
    syms eps_re_point real
    syms eps_im_point real

    gamma_vec = sym('gamma_n', [1 num_peaks]);
    sigma_vec = sym('sigma', [1 num_peaks]);
    w_0_vec = sym('w_0_', [1 num_peaks]);

    eps_re_str = "eps_infin";
    eps_im_str = "0";

    for func_num = 1:1:num_peaks
        eps_re_str = eps_re_str + " + (w_0_" + num2str(func_num) + "^2 * sigma" ...
            + num2str(func_num) + ") * (w_0_" + num2str(func_num) + "^2 - omega^2) / " ...
            + "((w_0_" + num2str(func_num) + "^2 - omega^2)^2 + omega^2 * " ...
            + "(gamma_n" + num2str(func_num) + " / (2 * pi))^2)";
        if(size(data,2)>2)
            eps_im_str = eps_im_str + " + (w_0_" + num2str(func_num) + "^2 * sigma" ...
            + num2str(func_num) + ") * (omega * (gamma_n" + num2str(func_num) + " / (2 * pi))) / " ...
            + "((w_0_" + num2str(func_num) + "^2 - omega^2)^2 + omega^2 * " ...
            + "(gamma_n" + num2str(func_num) + " / (2 * pi))^2)";
        end
    end
    eps_re_func = str2sym(eps_re_str);
    err_func_re = (eps_re_point - eps_re_func)^2;
    if (size(data,2)>2)
        eps_im_func = str2sym(eps_im_str);
        err_func_im = (eps_im_point - eps_im_func)^2;
        err_func_re = re_frac*err_func_re + im_frac*err_func_im;
    end
    
    % Finding Min Error
    % fprintf("Computing Gradient ...\n")
%     grad_err = sym(zeros((3*num_peaks+1),1));
%     grad_err(1) = diff(err_func_re,eps_infin);
    % fprintf("Index %d complete\n", 1)
    param_vec = sym(zeros((3*num_peaks+1),1));
    param_vec(1) = eps_infin;    
    for index = 1:1:num_peaks
        
%         grad_err(3*(index-1)+2) = diff(err_func_re,sigma_vec(index));
%         fprintf("Index %d complete\n", 3*(index-1)+2)
%         grad_err(3*(index-1)+3) = diff(err_func_re,w_0_vec(index));
%         fprintf("Index %d complete\n", 3*(index-1)+3)
%         grad_err(3*(index-1)+4) = diff(err_func_re,gamma_vec(index));
%         fprintf("Index %d complete\n", 3*(index-1)+4)

        param_vec(3*(index-1)+2) = sigma_vec(index);
        param_vec(3*(index-1)+3) = w_0_vec(index);
        param_vec(3*(index-1)+4) = gamma_vec(index);


    end
    param_guess = 10 * rand((3*num_peaks + 1),1); % Random initial values
    if (size(data,2)>2)
        err_func_num = sum(subs(err_func_re, {eps_re_point, omega, eps_im_point}, {eps_re_dat, freq_meep_dat, eps_im_dat}));
    else
        err_func_num = sum(subs(err_func_re, {eps_re_point, omega}, {eps_re_dat, freq_meep_dat}));
    end
    err_func_num = matlabFunction(err_func_num, 'Vars', {param_vec});
    lb = zeros((3*num_peaks + 1),1); % Lower bound (all parameters >= 0)
    ub = []; % No upper bound
    options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 1e6);
    fprintf("Finding Minima ...\n")
    warning('off', 'all')
    [best_param_opt, best_fval] = fmincon(err_func_num, param_guess, [], [], [], [], lb, ub, [], options);
%     for trial = 1:10
%         fprintf("Trial: %d\n", trial)
%         param_guess = 10 * rand((3*num_peaks + 1),1); % Random initial values
%         gs = GlobalSearch;
%         problem = createOptimProblem('fmincon', 'objective', err_func_num, ...
%         'x0', param_guess, 'lb', lb, 'ub', ub, 'options', options);
%         [param_opt_temp, fval_temp] = run(gs, problem);
% 
%         if fval_temp < best_fval
%             fprintf("New Minimum: %.2e\n", fval_temp)
%             best_fval = fval_temp;
%             best_param_opt = param_opt_temp;
%         end
%     end
    for trial = 1:500
        fprintf("Trial: %d\n", trial)
        param_guess = 20 * rand((3*num_peaks + 1),1); % Random starting values
        [param_opt_temp, fval_temp] = fmincon(err_func_num, param_guess, [], [], [], [], lb, ub, [], options);

        if fval_temp < best_fval
            fprintf("New Minimum: %.2e\n", fval_temp)
            best_fval = fval_temp;
            best_param_opt = param_opt_temp;
        end
    end
    warning('on', 'all')
    param_opt = best_param_opt;
    % fval = best_fval;
    View_Dispersion(param_opt(3:3:(3*num_peaks+1)), param_opt(4:3:(3*num_peaks+1)), param_opt(2:3:(3*num_peaks+1)),  ...
        param_opt(1), min(1./freq_meep_dat), max(1./freq_meep_dat), 1000, "BTO")
    figure(7)
    scatter((1./freq_meep_dat),eps_re_dat)
    title("Real Eps Data")
    xlabel("\lambda (um)")
    ylabel("Re Eps")
    if size(data,2)>2
        figure(8)
        scatter((1./freq_meep_dat),eps_im_dat)
        title("Imaginary Eps Data")
        xlabel("\lambda (um)")
        ylabel("Im Eps")
    end
    % crit_points = vpasolve(grad_err == 0, param_vec,rand*ones(10,1));
%     fprintf("Computing Hessian ...\n")
%     H = hessian(err_func,param_vec);
%     H_at_crit_txt = "subs(H, params_vec, ";
%     crit_points_txt = "[crit_points.eps_infin,";
%     for index = 1:1:num_peaks
%         crit_points_txt = crit_points_txt + " crit_points.sigma" ...
%             + num2str(index) + ", crit_points.w_0_" + num2str(index) ...
%             + ", crit_points.gamma_n" + num2str(index) + ",";
%     end
%     H_at_crit_txt = H_at_crit_txt + crit_points_txt + "])";
%     fprintf("Finding Hessian at Crit Points ...\n")
%     H_at_crit = str2sym(H_at_crit_txt);
%     fprintf("Solving for Eigen Values ...\n")
%     eig_vals = eig(H_at_crit);
% 
%     f_crit_txt = "subs(err_func, params_vec, " + crit_points_txt + "])";
%     f_crit = str2sym(f_crit_txt);
end