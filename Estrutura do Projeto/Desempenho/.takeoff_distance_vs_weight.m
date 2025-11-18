%% takeoff_distance_vs_weight_slope.m
% Distância de Decolagem × Peso com Efeito de Inclinação de Pista
% Base: Anderson / Roskam / Raymer
% Métodos: (1) Energy method (numérico), (2) Analytical approx

clear; clc; close all;
fprintf("=== TAKEOFF DISTANCE VS WEIGHT (with RUNWAY SLOPE) ===\n\n");

%% =================== INPUTS ====================
fprintf('Escolha modo de entrada:\n1 - Manual\n2 - Exemplo\n');
mode = input('Opção: ');

if mode == 1
    S       = input('Área da asa S [m^2]: ');
    AR      = input('Alongamento AR: ');
    e       = input('Eficiência Oswald e: ');
    CLmax_TO = input('CLmax em decolagem: ');
    CD0     = input('CD0: ');
    n_eng   = input('Nº motores: ');
    T_per_eng = input('Empuxo por motor ao nível do mar [N]: ');
    runway_available = input('Comprimento de pista disponível [m]: ');
    alt_analysis = input('Altitude da pista [m]: ');
    headwind = input('Vento de proa [m/s] (0 se nenhum): ');
    mu_r = input('Coeficiente de atrito de rolamento (0.02–0.04): ');
    slope_percent = input('Inclinação da pista [%]: (positivo = subida): ');
    W_min = input('Peso mínimo [kg]: ');
    W_max = input('Peso máximo [kg]: ');
    nW = input('Nº pontos de peso: ');
    h_obst = input('Altura do obstáculo [m] (ex: 15): ');
    fprintf('\nMétodo:\n1 - Energy (numérico)\n2 - Analytical (aprox)\n');
    method = input('Escolha método: ');
else
    % --------- Example default parameters ---------
    S = 16.2;
    AR = 8.5;
    e = 0.85;
    CLmax_TO = 2.0;
    CD0 = 0.028;
    n_eng = 1;
    T_per_eng = 9000;
    runway_available = 1500;
    alt_analysis = 0;
    headwind = 5;
    mu_r = 0.03;
    slope_percent = +2; % inclinação da pista (2% uphill)
    W_min = 2000; W_max = 4000; nW = 15;
    h_obst = 15;
    method = 1;
    fprintf('Usando exemplo completo com slope = +2%% (subida).\n');
end

%% --------- Slope conversion ----------
% slope_percent (%) → θ = arctan(slope/100)
theta = atan(slope_percent/100);

%% --------- Constants, ISA, thrust ----------
g = 9.80665;

atm = isa_calc(alt_analysis);
rho = atm.rho;

T_total_SL = n_eng * T_per_eng;

%% --------- Weight vector ----------
W_vec_kg = linspace(W_min, W_max, nW);

%% --------- Allocate ----------
S_ground = zeros(size(W_vec_kg));
S_climb  = zeros(size(W_vec_kg));
S_total  = zeros(size(W_vec_kg));
V_stall_vec = zeros(size(W_vec_kg));
V_rot_vec   = zeros(size(W_vec_kg));

%% ============================================================
%               MAIN LOOP — WEIGHT SWEEP
%% ============================================================
for i = 1:length(W_vec_kg)

    Wkg = W_vec_kg(i);
    WN  = Wkg*g;   % weight in N
    m   = Wkg;

    % Stall and rotation speeds
    V_stall = sqrt( (2*WN)/(rho*S*CLmax_TO) );
    V_rot   = 1.2 * V_stall;

    V_stall_vec(i) = V_stall;
    V_rot_vec(i)   = V_rot;

    % -------------------------------
    % METHOD 1 — ENERGY (NUMERICAL)
    % -------------------------------
    if method == 1
        V_air = linspace(0.1, V_rot, 500);
        q = 0.5*rho.*V_air.^2;

        % CL during ground roll
        CL_TO = 0.9 * CLmax_TO;
        L = q*S*CL_TO;

        CDi = (CL_TO.^2)./(pi*e*AR);
        CD = CD0 + CDi;
        D = q.*S.*CD;

        % constant thrust approx
        T = T_total_SL;

        % rolling friction
        F_roll = mu_r.*(WN - L);

        % ** SLOPE FORCE **
        F_slope = WN*sin(theta);

        % acceleration
        a = (T - D - F_roll)/m - g*sin(theta);
        a(a < 0.1) = 0.1;

        % ground speed (wind correction)
        V_ground = max(V_air - headwind, 0.1);

        % integrate
        s_arr = cumtrapz(V_air, V_ground./a);
        s_ground = s_arr(end);

        % climb to obstacle
        q_rot = 0.5*rho*V_rot^2;
        CL_req = WN/(q_rot*S);
        CD_rot = CD0 + (CL_req^2)/(pi*e*AR);
        D_rot = CD_rot*q_rot*S;

        grad = (T - D_rot)/WN;
        if grad <= 0
            s_climb = Inf;
        else
            ROC = grad*V_rot;
            t_climb = h_obst/ROC;
            s_climb = V_rot*t_climb;
        end

        S_ground(i) = s_ground;
        S_climb(i) = s_climb;
        S_total(i) = s_ground + s_climb;

    else
    % -------------------------------
    % METHOD 2 — Analytical (approx)
    % -------------------------------
        V_half = V_rot/2;
        q_half = 0.5*rho*V_half^2;
        CL_TO = 0.9*CLmax_TO;

        L_avg = q_half*S*CL_TO;
        CDi_avg = (CL_TO^2)/(pi*e*AR);
        D_avg = q_half*S*(CD0 + CDi_avg);

        F_roll_avg = mu_r*(WN - L_avg);

        % slope term
        a_avg = (T_total_SL - D_avg - F_roll_avg)/m - g*sin(theta);
        if a_avg < 0.1, a_avg = 0.1; end

        s_ground = V_rot^2 / (2*a_avg);

        % climb
        q_rot = 0.5*rho*V_rot^2;
        CL_req = WN/(q_rot*S);
        CD_rot = CD0 + (CL_req^2)/(pi*e*AR);
        D_rot = CD_rot*q_rot*S;
        grad = (T_total_SL - D_rot)/WN;

        if grad <= 0
            s_climb = Inf;
        else
            ROC = grad*V_rot;
            t_climb = h_obst/ROC;
            s_climb = V_rot*t_climb;
        end

        S_ground(i) = s_ground;
        S_climb(i) = s_climb;
        S_total(i) = s_ground + s_climb;
    end
end

%% ============================================================
%                           PLOTS
%% ============================================================
figure;
plot(W_vec_kg,S_ground,'-o','LineWidth',1.7); hold on;
plot(W_vec_kg,S_climb,'-s','LineWidth',1.7);
plot(W_vec_kg,S_total,'-d','LineWidth',2);
yline(runway_available,'--k','Runway');
xlabel('Peso [kg]'); ylabel('Distância [m]');
title(sprintf('Distância de Decolagem × Peso (Slope = %.1f%%)',slope_percent));
grid on;
legend('Ground roll','Climb','Total','Runway');

figure;
plot(W_vec_kg,V_stall_vec,'-o','LineWidth',1.6); hold on;
plot(W_vec_kg,V_rot_vec,'-s','LineWidth',1.6);
xlabel('Peso [kg]'); ylabel('Velocidade [m/s]');
title('V_{stall} e V_{rot} × Peso');
grid on; legend('V_{stall}','V_{rot}');

%% ---------- Output table ----------
T = table(W_vec_kg', V_stall_vec', V_rot_vec', S_ground', S_climb', S_total', ...
    'VariableNames',{'W_kg','V_stall','V_rot','S_ground','S_climb','S_total'});
disp(T);

%% ---------- Max weight for runway ----------
idx = find(S_total <= runway_available,1,'last');
if isempty(idx)
    fprintf("\nNenhum peso atende a pista disponível.\n");
else
    fprintf("\nPeso máximo dentro da pista disponível: %.1f kg  (S_total = %.1f m)\n", ...
        W_vec_kg(idx), S_total(idx));
end

%% ============================================================
% ISA Helper
%% ============================================================
function atm = isa_calc(h)

    g = 9.80665;
    R = 287.05;

    % Prealocação
    T = zeros(size(h));
    p = zeros(size(h));

    % Troposfera
    i1 = h <= 11000;
    T(i1) = 288.15 - 0.0065*h(i1);
    p(i1) = 101325 .* (T(i1)/288.15).^(g/(0.0065*R));

    % Estratosfera
    i2 = h > 11000;
    T(i2) = 216.65;
    p(i2) = 22632.06 .* exp(-g.*(h(i2)-11000)./(R*T(i2)));

    atm.T   = T;
    atm.p   = p;
    atm.rho = p ./ (R.*T);
    atm.a   = sqrt(1.4 .* R .* T);
end
