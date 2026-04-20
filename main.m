function dc_dc_mpc_fpga_fixed_point()
clc; clear; close all

%% SYSTEM PARAMETERS
Vin  = 24; Vref = 12; Rload = 10;
fs = 100e3; Ts = 1/fs;
L = 220e-6; C = 470e-6;

%% FPGA PARAMETERS
PWM_bits = 10; PWM_levels = 2^PWM_bits;

%% ADC PARAMETERS
ADC_bits = 12; ADC_levels = 2^ADC_bits;

%% SIMULATION
Tsim = 0.05; dt = Ts/200;
t = 0:dt:Tsim; N = length(t);

%% INIT VARIABLES
IL = zeros(1,N); Vout = zeros(1,N);
duty = zeros(1,N); pwm_signal = zeros(1,N);
duty(1) = 0.5;

%% LOSSES
Rds_on = 0.05; RL = 0.02;

%% MPC PARAMETERS
lambda = 0.05; duty_candidates = linspace(0.3,0.7,25);

%% FIXED-POINT PARAMETERS
duty_q_bits = PWM_bits; duty_q_levels = 2^duty_q_bits;

%% MAIN LOOP (FIXED-POINT MPC)
for k=1:N-1
    % Load step
    if t(k) > 0.025, Rload = 5; end

    % ADC quantization
    Vout_adc_q = round((Vout(k)/Vin)*ADC_levels)/ADC_levels*Vin;
    IL_adc_q   = round((IL(k)/10)*ADC_levels)/ADC_levels*10;

    % MPC Fixed-Point
    best_cost = 1e9; best_duty = duty(k);
    duty_candidates_q = round(duty_candidates*duty_q_levels)/duty_q_levels;

    for d = duty_candidates_q
        IL_pred = IL_adc_q + ((Vin*d - Vout_adc_q)/L)*Ts;
        V_pred  = Vout_adc_q + ((IL_adc_q - Vout_adc_q/Rload)/C)*Ts;
        cost = (Vref - V_pred)^2 + lambda*(d - duty(k))^2;
        if cost < best_cost
            best_cost = cost; best_duty = d;
        end
    end
    duty(k+1) = best_duty;

    % PWM signal generation
    duty_count = round(duty(k)*duty_q_levels);
    counter = mod(k,duty_q_levels);
    pwm_signal(k) = counter < duty_count;

    % DC-DC Converter Equations
    Vsw = pwm_signal(k)*Vin;
    dIL = (Vsw - Vout(k) - IL(k)*RL)/L;
    dVout = (IL(k) - Vout(k)/Rload)/C;
    IL(k+1) = IL(k) + dIL*dt;
    Vout(k+1) = Vout(k) + dVout*dt;
end

%% PERFORMANCE
steady = round(N*0.8):N;
Iout = Vout./Rload; Iin = duty.*IL;
Pin  = Vin*mean(Iin(steady));
Pout = mean(Vout(steady).*Iout(steady));
Pmosfet = mean((IL(steady).^2)*Rds_on);
PL = mean((IL(steady).^2)*RL);
Ploss = Pmosfet + PL;
efficiency = (Pout/(Pout+Ploss))*100;
Vout_avg = mean(Vout(steady));
ripple = max(Vout(steady)) - min(Vout(steady));
settling_index = find(abs(Vout-Vref)<0.02*Vref,1);
settling_time = t(settling_index);
Pin_inst  = Vin.*Iin; Pout_inst = Vout.*Iout;

%% PRINT FIXED-POINT OPTIMIZATION RESULTS
fprintf('\n==== FIXED-POINT ARITHMETIC OPTIMIZATION RESULTS ====\n');
fprintf('PWM Bits: %d, ADC Bits: %d\n', PWM_bits, ADC_bits);
fprintf('Duty Cycle Range: %.4f - %.4f\n', min(duty), max(duty));
fprintf('Efficiency: %.2f %%\n', efficiency);
fprintf('Voltage Ripple: %.4f V\n', ripple);
fprintf('Settling Time: %.6f s\n', settling_time);

%% ===========================================
% FIXED-POINT OPTIMIZATION PLOTS
%% ===========================================

% 1. Duty Cycle Evolution
figure
plot(t,duty,'b','LineWidth',2)
title('Fixed-Point MPC: Duty Cycle Evolution')
xlabel('Time (s)'); ylabel('Duty Cycle'); grid on

% 2. Output Voltage vs Reference
figure
plot(t,Vout,'b','LineWidth',2)
hold on
plot(t,Vref*ones(size(t)),'r--','LineWidth',2)
title('Fixed-Point MPC: Output Voltage Tracking')
xlabel('Time (s)'); ylabel('Voltage (V)')
legend('Vout','Vref'); grid on

% 3. Steady-State Ripple
figure
plot(t(steady),Vout(steady),'LineWidth',2)
title(['Fixed-Point MPC: Voltage Ripple = ',num2str(ripple*1000,'%.2f'),' mV'])
xlabel('Time (s)'); ylabel('Voltage (V)'); grid on

% 4. Efficiency vs Load
load_values = [20 15 10 7 5 3]; eff_values = zeros(size(load_values));
for i=1:length(load_values)
    Rtest = load_values(i);
    Iout_test = Vout ./ Rtest;
    Pout_test = mean(Vout(steady).*Iout_test(steady));
    eff_values(i) = (Pout_test/(Pout_test+Ploss))*100;
end
figure
plot(load_values, eff_values,'-o','LineWidth',2)
title('Fixed-Point MPC: Efficiency vs Load')
xlabel('Load Resistance (Ohms)'); ylabel('Efficiency (%)'); grid on

% 5. Power Loss Analysis
Pmosfet_inst = (IL.^2)*Rds_on;
PL_inst      = (IL.^2)*RL;
TotalLoss    = Pmosfet_inst + PL_inst;
figure
plot(t,Pmosfet_inst,'LineWidth',1.5)
hold on
plot(t,PL_inst,'LineWidth',1.5)
plot(t,TotalLoss,'LineWidth',2)
title('Fixed-Point MPC: Power Loss Breakdown')
xlabel('Time (s)'); ylabel('Power Loss (W)')
legend('MOSFET Loss','Inductor Loss','Total Loss'); grid on

%% ===============================
% MPC Cost Convergence Tracking
%% ===============================
cost_history = zeros(1,N);  % store best cost at each timestep

for k=1:N-1
    % Load change
    if t(k) > 0.025
        Rload = 5;
    end

    % ADC quantization
    Vout_adc_q = round((Vout(k)/Vin)*ADC_levels)/ADC_levels*Vin;
    IL_adc_q   = round((IL(k)/10)*ADC_levels)/ADC_levels*10;

    % MPC Fixed-Point
    best_cost = 1e9; best_duty = duty(k);
    duty_candidates_q = round(duty_candidates*duty_q_levels)/duty_q_levels;

    for d = duty_candidates_q
        IL_pred = IL_adc_q + ((Vin*d - Vout_adc_q)/L)*Ts;
        V_pred  = Vout_adc_q + ((IL_adc_q - Vout_adc_q/Rload)/C)*Ts;
        cost = (Vref - V_pred)^2 + lambda*(d - duty(k))^2;

        if cost < best_cost
            best_cost = cost;
            best_duty = d;
        end
    end

    duty(k+1) = best_duty;
    cost_history(k) = best_cost;  % save best cost for convergence plot

    % PWM generation and converter
    duty_count = round(duty(k)*duty_q_levels);
    counter = mod(k,duty_q_levels);
    pwm_signal(k) = counter < duty_count;

    Vsw = pwm_signal(k)*Vin;
    dIL = (Vsw - Vout(k) - IL(k)*RL)/L;
    dVout = (IL(k) - Vout(k)/Rload)/C;
    IL(k+1) = IL(k) + dIL*dt;
    Vout(k+1) = Vout(k) + dVout*dt;
end

%% ===============================
% Convergence Plot
%% ===============================
figure
plot(t, cost_history, 'LineWidth',2)
title('MPC Cost Function Convergence (Fixed-Point)')
xlabel('Time (s)')
ylabel('MPC Cost')
grid on