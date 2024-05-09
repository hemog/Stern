% MATLAB kodu
% Gerekli kütüphanelerin eklenmesi
addpath('pyXSteam'); % pyXSteam kütüphanesinin MATLAB sürümü eklenmeli

% Constants
P6 = 0.5;
T7_values = 448:50:873;

% Fonksiyonlar
function [eta_th, m, w_net, q_net, w_p, w_t, q_L, q_H] = calculate_eta_th(P6, T7_value)
    steam_table = XSteam(XSteam.UNIT_SYSTEM_BARE);
    
    state_point_1.P = 0.0075; % MPa
    state_point_1.h = steam_table.hL_p(state_point_1.P);
    state_point_1.s = steam_table.sL_p(state_point_1.P);
    state_point_1.T = steam_table.tsat_p(state_point_1.P);

    state_point_2.P = P6;
    state_point_2.T = state_point_1.T;
    state_point_2.h = steam_table.hL_t(state_point_2.T);
    state_point_2.s = steam_table.sL_t(state_point_2.T);

    state_point_3.P = state_point_2.P;
    state_point_3.h = steam_table.hL_p(state_point_3.P);
    state_point_3.s = steam_table.sL_p(state_point_3.P);

    state_point_4.P = 9.5;
    state_point_4.s = state_point_3.s;
    state_point_4.T = steam_table.t_ps(state_point_4.P, state_point_4.s);
    state_point_4.h = steam_table.hL_t(state_point_4.T);

    state_point_5.P = state_point_4.P;
    state_point_5.T = 873;
    state_point_5.h = steam_table.h_pt(state_point_5.P, state_point_5.T);
    state_point_5.s = steam_table.s_pt(state_point_5.P, state_point_5.T);

    eta_hpt = 1;
    state_point_6.P = state_point_3.P;
    state_point_6.s_s = state_point_5.s;
    state_point_6.h_s = steam_table.h_ps(state_point_6.P, state_point_6.s_s);
    state_point_6.h = state_point_5.h - eta_hpt * (state_point_5.h - state_point_6.h_s);
    state_point_6.x = steam_table.x_ph(state_point_6.P, state_point_6.h);
    state_point_6.T = steam_table.t_ps(state_point_6.P, state_point_5.s);

    state_point_7.P = state_point_6.P;
    state_point_7.T = T7_value;
    state_point_7.h = steam_table.h_pt(state_point_7.P, state_point_7.T);
    state_point_7.s = steam_table.s_pt(state_point_7.P, state_point_7.T);

    state_point_8.P = state_point_1.P;
    eta_lpt = 1;
    state_point_8.s = state_point_7.s;
    state_point_8.h = steam_table.h_ps(state_point_8.P, state_point_8.s);
    state_point_8.x = steam_table.x_ps(state_point_8.P, state_point_8.s);

    m = ((state_point_3.h - state_point_2.h)) / ((state_point_6.h) - state_point_2.h);

    w_p1 = (state_point_2.h - state_point_1.h) * (1 - m);
    w_p2 = (state_point_4.h - state_point_3.h);
    w_p = w_p1 + w_p2;

    w_hpt = (state_point_5.h - state_point_6.h);
    w_lpt = (state_point_7.h - state_point_8.h) * (1 - m);
    w_t = w_hpt + w_lpt;

    w_net = (w_hpt + w_lpt) - (w_p1 + w_p2);

    q_H = state_point_5.h - state_point_4.h + (state_point_7.h - state_point_6.h) * (1 - m);
    q_L = (state_point_8.h - state_point_1.h) * (1 - m);
    q_net = q_H - q_L;

    eta_th = 1 - (q_L / q_H);
end

function fval = objective_function(params)
    P6 = params(1);
    T7 = params(2);
    [eta_th, ~, ~, ~, ~, ~, ~, ~] = calculate_eta_th(P6, T7);
    fval = -eta_th;
end

% Optimize
initial_guess = [0.5, 673];
lb = [0.0075, 448];
ub = [9.5, 873];
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1000);
[result, max_eta_th] = fmincon(@objective_function, initial_guess, [], [], [], [], lb, ub, [], options);
optimal_P6 = result(1);
optimal_T7 = result(2);

% Sonuçları göster
fprintf('Optimized T7: %.2f K\n', optimal_P6);
fprintf('Optimized P6: %.2f MPa\n', optimal_T7);
fprintf('Maximum eta_th: %.4f\n', -max_eta_th);
