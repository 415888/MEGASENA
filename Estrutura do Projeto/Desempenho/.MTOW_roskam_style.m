%% estimate_MTOW_roskam_style.m
% Estima MTOW via método iterativo (estilo Roskam/Raymer) com inputs.
% Resultado: MTOW que atende requisitos de missão e restrições.
% Autor: ChatGPT - versão adaptada (use como base; ajuste coeficientes para seu projeto)

clear; clc; close all;
fprintf('=== MTOW Estimation (Roskam-style iterative) ===\n\n');

%% -------------------------
% USER INPUTS (mission / geometry / propulsion / requirements)
%% -------------------------
% Geometry / weights
S = input('Área da asa S [m^2] (ex: 16.2): ');
AR = input('Alongamento AR (ex: 8.5): ');
e = input('Eficiência Oswald e (ex: 0.85): ');
CLmax_TO = input('CL_max (takeoff, com flaps) (ex: 2.0): ');
CLmax_clean = input('CL_max clean (ex: 1.6): ');

payload_kg = input('Carga útil (payload) [kg] (ex: 500): ');
crew_weight_kg = input('Peso tripulação [kg] (total) (ex: 100): ');

% Mission
Range_km = input('Alcance de missão (único trecho) [km] (ex: 1000): ');
Reserve_fraction = input('Reserva combustível (fraction of trip fuel) (ex: 0.05): ');
taxi_fraction = input('Taxi fuel fraction of takeoff (ex: 0.01): ');
climb_fraction = input('Climb fraction (ex: 0.02): ');

% Engine / propulsion
engine_type = input('Engine type: 1=Turbojet/TP, 2=Piston/Eletrico (enter 1 or 2): ');
if engine_type==1
    % turboprop/turbofan use specific fuel consumption in N/W or 1/s?
    sfc = input('SFC (1/hr) for cruise (ex: 0.5 for turboprop, units 1/hr): '); % 1/hr typical
    sfc = sfc/3600; % convert to 1/s
else
    % piston/electric: user can input fuel energy or battery energy assumptions
    sfc = input('Equivalent "sfc" (1/hr) or energy-based parameter (ex: 0.3) : ');
    sfc = sfc/3600;
end
eta_prop = input('Eficiência propulsora (0-1) ex: 0.85: ');

% Power / thrust sizing guess
T_over_W_guess = input('Chute inicial T/W (thrust per weight) (ex: 0.3 N per N = 0.3): ');
P_over_W_guess = input('Chute P/W [W/N] (potencia disponivel por Newton de peso), se aplicável (ex: 0.08 W/N): ');
% runway and climb constraints
runway_length = input('Pista disponível [m] (ex: 1500): ');
h_obstacle = input('Altura do obstáculo [m] (ex: 15): ');

% Empty weight correlation coefficients (Roskam-style empirical)
fprintf('\nEmpty weight estimation: W_empty = A * MTOW^B\n');
A = input('Coeficiente A (ex: 0.5): ');
B = input('Expoente B (ex: 0.9): ');

% Structural / regulatory
V_stall_limit = input('V_stall limit [m/s] (ex: 35) - leave empty for no explicit limit: ');
if isempty(V_stall_limit), V_stall_limit = Inf; end

% Iteration control
MTOW_initial = input('Chute inicial MTOW [kg] (ex: 6000): ');
if isempty(MTOW_initial), MTOW_initial = 6000; end
max_iter = input('Max iterations (ex: 50): ');
if isempty(max_iter), max_iter = 40; end
tol = input('Tolerância convergência MTOW (kg) (ex: 1.0): ');
if isempty(tol), tol = 1.0; end

%% -------------------------
% constants and defaults
g = 9.81;
rho0 = 1.225;
a0 = 340; % m/s approx speed of sound at sea level
% Derived (aerodynamic)
CL_alpha = 2*pi*0.9; % rough 2D to 3D correction approx; unused directly
% mission conversions
Range = Range_km * 1000;

% set limits and margin factors
safety_margin_takeoff = 1.15; % require available runway <= runway_length / margin
safety_margin_climb = 0.98;  % require climb gradient margin
% Sizing parameters
MTOW = MTOW_initial; % kg initial
fprintf('\nStarting iteration for MTOW (initial %.1f kg)\n', MTOW);

%% -------------------------
% iterative loop: guess MTOW, compute components, check constraints, update
%% -------------------------
converged = false;
history = [];

for it = 1:max_iter
    % compute weights
    MTOW_N = MTOW * g;
    % empty weight via empirical correlation
    W_empty = A * (MTOW^B);   % kg
    W_payload = payload_kg + crew_weight_kg; % kg (simple)
    
    % fuel estimation: iterative because fuel depends on MTOW (weight changes)
    % Use Breguet for cruise fuel fraction: W_final = W_initial * exp(-Range * c / (V * L/D))
    % We need assumed cruise speed V_cruise and L/D cruise. We'll ask user (with defaults)
    if it==1
        V_cruise = input('Velocidade de cruzeiro assumida [m/s] (ex: 200): ');
        if isempty(V_cruise), V_cruise = 200; end
        LD_cruise = input('L/D médio de cruzeiro (ex: 12): ');
        if isempty(LD_cruise), LD_cruise = 12; end
    end
    
    % approximate c = sfc * g? For Breguet with power specific fuel consumption in 1/s:
    % Breguet (mass-based): Wf/Wi = 1 - exp(-Range * c / (V*L/D))
    % Here c should be fuel consumption rate per unit weight (1/s). Use sfc (1/s) * g?
    % We'll compute c_eff = sfc * g (approx) so units consistent: sfc [1/s] * g [N/kg] => 1/s? 
    % To keep practical: use c_eff = sfc (1/s) as user input; may need tuning.
    c_eff = sfc; % user supplied earlier in 1/s
    if c_eff <= 0
        error('SFC (c_eff) deve ser > 0');
    end
    
    % compute cruise fuel fraction (mass fraction consumed during cruise)
    % W_after_cruise = W_before_cruise * exp(-Range * c_eff / (V_cruise * LD_cruise));
    frac_cruise = 1 - exp(-Range * c_eff / (V_cruise * LD_cruise)); % fraction consumed during cruise
    if frac_cruise < 0, frac_cruise = 0; end
    % Other mission fractions (taxi, takeoff, climb, landing). Use simple sums (fractions of takeoff weight).
    taxi_frac = taxi_fraction; % fraction of takeoff weight
    climb_frac = climb_fraction;
    reserve_frac = Reserve_fraction;
    % total fuel fraction approximately:
    % Wf_total = MTOW * (taxi_frac + climb_frac + frac_cruise + reserve_frac)
    fuel_mass_kg = MTOW * (taxi_frac + climb_frac + frac_cruise + reserve_frac);
    % ensure fuel mass non-negative
    if fuel_mass_kg < 0, fuel_mass_kg = 0; end
    
    W_fuel = fuel_mass_kg; % kg
    W_empty_mass = W_empty; % kg
    W_payload_mass = W_payload; % kg
    W_cg_payload = W_payload_mass; % used later
    
    % zero-fuel weight (ZFW)
    W_zfw = W_empty_mass + W_payload_mass;
    % takeoff gross weight check
    W_takeoff = W_zfw + W_fuel;
    
    % compute T/W and P/W based on guesses (user inputs). We can update T based on MTOW.
    % If user provided T/W guess, compute total thrust available or power available:
    T_over_W = T_over_W_guess; % ratio (dimensionless)
    T_total_N = T_over_W * (MTOW_N); % N total thrust
    % Power available approx
    P_over_W = P_over_W_guess; % W per N
    P_avail_total = P_over_W * MTOW_N; % W
    
    % Aerodynamic derived
    % Stall speed
    rho = rho0; % sea-level baseline (we can compute at altitude later)
    V_stall = sqrt((2 * (MTOW * g)) / (rho * S * CLmax_clean)); % m/s
    
    % Takeoff performance approximate:
    % We'll compute ground roll using energy method with average acceleration:
    % a_avg = (T - D_avg - mu*(W - L_avg)) / m
    mu_roll = 0.03; % rolling friction
    m = MTOW; % kg
    % approximate D_avg and L_avg during roll using CL_TO (fraction of CLmax_TO)
    CL_TO = 0.9 * CLmax_TO;
    V_rot = 1.2 * sqrt((2 * MTOW * g) / (rho * S * CLmax_TO));
    V_array = linspace(0, V_rot, 200);
    q_roll = 0.5 * rho .* V_array.^2;
    L_roll = CL_TO .* q_roll .* S;
    CD_TO = CD0_from_AR(AR, e) + (CL_TO.^2)/(pi*e*AR); % CD0 estimate from AR (helper creates CD0 if needed)
    D_roll = CD_TO .* q_roll .* S;
    T_roll = T_total_N * ones(size(V_array));
    a_roll = (T_roll - D_roll - mu_roll*(MTOW*g - L_roll))/m; % m/s^2
    a_roll(a_roll<0.1) = 0.1; % floor
    s_roll = cumtrapz(V_array ./ a_roll);
    ground_run = s_roll(end);
    % airborne distance to obstacle (approx)
    % use climb rate at V_rot: ROC_at_vrot ~ (T - D)/W * V ~ (P_av - P_req)/W ; use simple approach
    % compute D at V_rot
    q_vrot = 0.5 * rho * V_rot^2;
    CL_req_vrot = (MTOW*g) / (q_vrot * S);
    CD_vrot = CD0_from_AR(AR, e) + (CL_req_vrot^2)/(pi*e*AR);
    D_vrot = CD_vrot * q_vrot * S;
    % thrust at Vrot
    T_at_vrot = T_total_N;
    % approximate climb gradient
    grad_climb = (T_at_vrot - D_vrot)/(MTOW*g); % dimensionless
    if grad_climb <= 1e-6
        climb_dist = Inf;
    else
        ROC_vrot = grad_climb * V_rot; % m/s
        t_to_clear = h_obstacle / ROC_vrot;
        climb_dist = V_rot * t_to_clear;
    end
    S_to_total = ground_run + climb_dist;
    
    % Stall constraint:
    meets_stall = V_stall <= V_stall_limit;
    
    % Climb gradient constraint per FAR-25/AEO typical:
    % require gradient >= 0.024 (2.4%) for multi-engine AEO at V2 approx; user can adjust
    required_grad = 0.024; % default AEO threshold
    meets_climb = grad_climb >= required_grad;
    
    % Runway constraint (with margin)
    meets_runway = (S_to_total * safety_margin_takeoff) <= runway_length;
    
    % Structural check simplified: wing loading vs CLmax gives maneuvering factor; here we just check L/W * n_limit
    wing_loading = (MTOW*g)/S; % N/m2
    % simplistic g-limit: ensure wings produce required lift at limit g (if known). We'll skip complex structural calc.
    meets_struct = true;
    
    % collect flags
    ok_flags = [meets_runway, meets_climb, meets_stall, meets_struct];
    
    % Logging history
    history = [history; it, MTOW, W_empty, W_fuel, S_to_total, grad_climb, V_stall, meets_runway, meets_climb];
    
    % Decide adjustment direction:
    if all(ok_flags)
        % candidate MTOW ok -> try increasing to see margin (or converge)
        % but our goal is the minimum MTOW that satisfies constraints for required payload+fuel.
        % So we seek fixed point where MTOW = W_empty(MTOW) + payload + fuel(MTOW)
        % compute implied MTOW_from_components
        MTOW_implied = W_empty + W_payload_mass + W_fuel;
        diff = MTOW_implied - MTOW;
        % if difference small -> converged
        fprintf('Iter %d: MTOW=%.1f kg OK. implied MTOW=%.1f kg diff=%.2f kg\n', it, MTOW, MTOW_implied, diff);
        if abs(diff) <= tol
            converged = true;
            MTOW_final = MTOW;
            break;
        else
            % Update MTOW toward implied (simple relaxation)
            MTOW = 0.6*MTOW + 0.4*MTOW_implied;
        end
    else
        % not ok: need to reduce MTOW (if runway/climb/stall fail) or increase (if fuel insufficient)
        MTOW_implied = W_empty + W_payload_mass + W_fuel;
        diff = MTOW_implied - MTOW;
        fprintf('Iter %d: MTOW=%.1f kg NOT OK. runway?%d climb?%d stall?%d | implied MTOW=%.1f diff=%.2f\n', ...
            it, MTOW, meets_runway, meets_climb, meets_stall, MTOW_implied, diff);
        % if runway/climb/stall fail -> reduce MTOW by factor
        if (~meets_runway) || (~meets_climb) || (~meets_stall)
            MTOW = MTOW * 0.95; % reduce 5%
        else
            % otherwise move toward implied
            MTOW = 0.6*MTOW + 0.4*MTOW_implied;
        end
    end
    
    % safety: avoid negative/too small MTOW
    if MTOW < 100, error('MTOW fell below 100 kg - check inputs'); end
end

% After iterations:
if ~exist('MTOW_final','var')
    MTOW_final = MTOW;
    warning('Convergence não alcançada estritamente dentro de iterações. Resultado aproximado.');
end

%% -------------------------
% OUTPUT SUMMARY & PLOTS
%% -------------------------
fprintf('\n=== RESULTADO FINAL ===\n');
fprintf('MTOW final aproximado: %.2f kg\n', MTOW_final);
fprintf('Empty weight (W_empty): %.2f kg\n', W_empty);
fprintf('Fuel mass (kg): %.2f kg\n', W_fuel);
fprintf('ZFW (kg): %.2f kg\n', W_zfw);
fprintf('Takeoff distance (m): %.1f (ground %.1f + climb %.1f)\n', S_to_total, ground_run, climb_dist);
fprintf('Climb gradient at V_rot: %.3f (%.2f%%)\n', grad_climb, 100*grad_climb);
fprintf('V_stall (m/s): %.2f\n', V_stall);

% plot history: MTOW vs iteration, S_to vs MTOW, grad vs MTOW
history_table = array2table(history,'VariableNames',{'iter','MTOW','W_empty','W_fuel','S_to','grad_climb','V_stall','runwayOK','climbOK'});
figure; subplot(2,2,1);
plot(history(:,1), history(:,2), '-o'); xlabel('Iter'); ylabel('MTOW [kg]'); title('MTOW vs Iter');

subplot(2,2,2);
plot(history(:,1), history(:,5), '-o'); xlabel('Iter'); ylabel('S_{TO} [m]'); title('Takeoff distance history');

subplot(2,2,3);
plot(history(:,1), history(:,6), '-o'); xlabel('Iter'); ylabel('grad climb'); title('Climb gradient history');

subplot(2,2,4);
plot(history(:,1), history(:,7), '-o'); xlabel('Iter'); ylabel('V_{stall} [m/s]'); title('Vstall history');

% S_to vs MTOW parametric sweep for visualization
MTOW_sweep = linspace(max(100,MTOW_final*0.7), MTOW_final*1.3, 30);
S_to_sweep = zeros(size(MTOW_sweep));
Vstall_sweep = zeros(size(MTOW_sweep));
grad_sweep = zeros(size(MTOW_sweep));
for i=1:length(MTOW_sweep)
    mt = MTOW_sweep(i);
    rho = rho0;
    V_st = sqrt((2 * (mt * g)) / (rho * S * CLmax_clean));
    Vrot = 1.2 * sqrt((2 * mt * g) / (rho * S * CLmax_TO));
    % approximate ground run using simple acceleration a = (T - avgD - mu*(W-L))/m
    CL_TO = 0.9*CLmax_TO;
    Varr = linspace(0, Vrot, 150);
    qv = 0.5*rho.*Varr.^2;
    L_v = CL_TO .* qv .* S;
    CDt = CD0_from_AR(AR,e) + (CL_TO.^2)/(pi*e*AR);
    D_v = CDt .* qv .* S;
    Tt = T_over_W_guess * (mt * g);
    aavg = mean((Tt - D_v - mu_roll*(mt*g - L_v))/mt);
    if aavg < 0.1, aavg = 0.1; end
    S_to_sweep(i) = Vrot^2/(2*aavg) + (Vrot*(h_obstacle / max(0.1,( (Tt - interp1(Varr, D_v, Vrot, 'linear', D_v(end)))/(mt*g) * Vrot))));
    Vstall_sweep(i) = V_st;
    grad_sweep(i) = (Tt - interp1(Varr, D_v, Vrot, 'linear', D_v(end)))/(mt*g);
end

figure;
yyaxis left; plot(MTOW_sweep, S_to_sweep,'-o'); ylabel('S_{TO} [m]');
yyaxis right; plot(MTOW_sweep, Vstall_sweep,'-x'); ylabel('V_{stall} [m/s]');
xlabel('MTOW [kg]'); title('Sensibilidade S_{TO} e V_{stall} ao MTOW'); grid on;

figure;
plot(MTOW_sweep, grad_sweep); xlabel('MTOW [kg]'); ylabel('Climb gradient'); grid on;
title('Gradiente de subida vs MTOW');

fprintf('\nHistórico final salvo em variable "history" no workspace.\n');

%% Helper function: estimate CD0 from AR & Re (simple empirical)
function CD0 = CD0_from_AR(AR,e)
    % Empirical baseline CD0 estimate based on aspect ratio (very rough)
    % For more accuracy, user should input CD0 directly.
    CD0 = 0.02 + 0.01*(6/AR); % simple scaling: AR high -> lower CD0
    % ensure positive
    if CD0 < 0.01, CD0 = 0.01; end
end
