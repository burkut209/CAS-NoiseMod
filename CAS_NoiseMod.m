%% ================================================================
%  CAS-NoiseMod (WCL2025-2845) - UNIFIED MATLAB SCRIPT
%  Author: Serdar Prencuva (clean unified version)
%
%  Purpose:
%   (A) CAS-NoiseMod simulation: Bob BER, Eve BER (mismatched threshold)
%   (B) Bob theoretical BEP via moment-matched Gamma averaging
%   (C) Conventional NoiseMod (no pre-scaling)
%       - Bob BER baseline
%       - Eve BER baseline
%   (D) Figures:
%       - Fig. 4(a): Bob BER vs delta for multiple alpha (simulation + theory)
%       - Fig. 4(b): 3D visualization of Bob and Eve BER vs delta and alpha
%       - Fig. 4(c): 2D comparison (Bob/Eve CAS vs Bob/Eve Baseline) at fixed alpha
%
%  Notes:
%   - CAS-NoiseMod model: s[n] = rho * x[n], rho = |h| (dual-channel effect),
%     y[n] = h*s[n] + w[n]. See System Model Sec. II.  (Eq. (3)-(6)) 
%   - Conventional NoiseMod baseline: s[n] = x[n], y[n] = h*x[n] + w[n]
%   - Eve in CAS-NoiseMod: y_e[n] = h_e*s[n] + w_e[n],
%     but Eve does not know rho = |h| and thus uses a mismatched threshold.
%   - Random seed fixed for repeatability.
% ================================================================

clear; close all; clc;

%% ================================================================
%  SECTION 1: GLOBAL PARAMETERS (paper-consistent knobs)
% ================================================================
rng(3, 'twister');                  % Fixed seed for repeatability

N        = 120;                     % Samples per bit (paper uses N=120)
P        = 1;                       % Average transmit power
nbits_F4 = 1e6;                     % Bits for Fig. 4(a) style curves
nbits_F5 = 1e6;                     % Bits for Fig. 4(c) comparison plot
nbits_3D = 5e5;                     % Bits for 3D surface plot (reduced for speed)

delta_dB_F4 = -10:1:40;             % delta(dB) sweep for Fig. 4(a)
delta_F4    = 10.^(delta_dB_F4/10);

alpha_list_F4 = [1, 1.5, 3, 10];    % Fig. 4(a) alpha set

% Fig. 4(c) parameters (2D comparison plot)
delta_dB_F5 = -10:1:45;             % delta(dB) sweep for Fig. 4(c)
delta_F5    = 10.^(delta_dB_F5/10);
alpha_F5    = 3.5;                  % Fixed alpha for Fig. 4(c)

% Fig. 4(b) parameters (3D surface plot)
delta_dB_3D = -10:5:40;             % Coarser delta sweep for 3D plot
alpha_3D    = 1:0.5:3;              % Alpha range for 3D plot

% Toggle baseline overlay for comparison
INCLUDE_BASELINE_NOISEMOD = true;

%% ================================================================
%  SECTION 2: FIG. 4(a) - Bob BER vs delta for multiple alpha
%             Simulation + Theoretical (CAS-NoiseMod / Bob only)
%             Shows BER performance across different variance ratios.
% ================================================================
BER_Bob_sim_F4  = zeros(numel(alpha_list_F4), numel(delta_F4));
BER_Bob_the_F4  = zeros(numel(alpha_list_F4), numel(delta_F4));

fprintf('=== Fig. 4(a): Bob BER (CAS-NoiseMod) for multiple alpha ===\n');

% Generate one bit stream per alpha (kept independent across alpha for clarity)
for ia = 1:numel(alpha_list_F4)
    alpha = alpha_list_F4(ia);

    % Variances (Eq. (2) style):
    % sigma_L^2 = 2P/(1+alpha), sigma_H^2 = alpha*sigma_L^2
    sigmaL2 = 2*P/(1+alpha);
    sigmaH2 = alpha*sigmaL2;

    bits = randi([0 1], 1, nbits_F4);

    fprintf('  alpha = %.2f: sim...\n', alpha);

    for iD = 1:numel(delta_F4)
        delta = delta_F4(iD);
        sigmaW2 = sigmaL2 / delta;     % matches definition delta = sigmaL2^2 / sigmaW2

        err = 0;

        for k = 1:nbits_F4
            b = bits(k);

            % Legitimate channel (Rayleigh): h ~ CN(0,1), rho=|h|
            h  = (randn + 1j*randn)/sqrt(2);
            rho = abs(h);

            % NoiseMod symbol: x[n] real Gaussian with variance sigma_b^2 (Eq. (1))
            if b == 0
                x = sqrt(sigmaL2) * randn(1, N);
            else
                x = sqrt(sigmaH2) * randn(1, N);
            end

            % CAS-NoiseMod pre-scaling: s[n] = rho * x[n] (Eq. (3))
            s = rho * x;

            % AWGN at Bob: w ~ CN(0, sigmaW2)
            w = sqrt(sigmaW2/2) * (randn(1, N) + 1j*randn(1, N));

            % Received: y = h*s + w (Eq. (5))
            y = h * s + w;

            % Decision statistic: sample variance estimate (Eq. (7))
            Shat = mean(abs(y).^2);

            % Threshold (equal-shape / variance-threshold form; Eq. (8)/(16) spirit)
            % Under CAS, Bob's conditional variances: Cb = (rho^4)*sigma_b^2 + sigmaW2
            C0 = (rho^4) * sigmaL2 + sigmaW2;
            C1 = (rho^4) * sigmaH2 + sigmaW2;
            gamma = (C1*C0)/(C1-C0) * log(C1/C0);

            bhat = (Shat > gamma);
            err  = err + xor(bhat, b);
        end

        BER_Bob_sim_F4(ia, iD) = err / nbits_F4;
    end

    fprintf('  alpha = %.2f: theory...\n', alpha);
    BER_Bob_the_F4(ia, :) = BEP_Bob_CAS_Theory(delta_dB_F4, N, alpha, P);
end

% Handle alpha=1 case: variances identical -> random guessing (BER=0.5)
if any(alpha_list_F4 == 1)
    idx1 = find(alpha_list_F4 == 1, 1);
    BER_Bob_the_F4(idx1, :) = 0.5;
end

%% Plot Fig. 4(a) - Bob BER vs delta for multiple alpha values
figure('Color','white','Position',[100 100 850 600]); hold on; grid on;
set(gca,'YScale','log');

% Plot simulations for each alpha with markers
mk = {'o','s','^','d'};
for ia = 1:numel(alpha_list_F4)
    semilogy(delta_dB_F4, BER_Bob_sim_F4(ia,:), ['-' mk{ia}], 'LineWidth', 1.6, 'MarkerSize', 4, ...
        'DisplayName', sprintf('\\alpha = %.1f (Sim)', alpha_list_F4(ia)));
end
% Plot one "Theoretical" style (black dashed) for each alpha
for ia = 1:numel(alpha_list_F4)
    semilogy(delta_dB_F4, BER_Bob_the_F4(ia,:), '--', 'LineWidth', 1.4, 'Color', 'k', ...
        'HandleVisibility','off');
end
% Add a dummy theoretical handle for legend cleanliness
hTheo = semilogy(nan, nan, '--k', 'LineWidth', 1.4, 'DisplayName', 'Theoretical');

xlabel('\delta (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('BER', 'FontSize', 14, 'FontName', 'Times New Roman');
ylim([1e-3 1]);
xlim([min(delta_dB_F4) max(delta_dB_F4)]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
set(gca, 'Box', 'on', 'LineWidth', 1.5);
legend('Location','northeast');
title('(a)', 'FontSize', 14, 'FontName', 'Times New Roman');

print(gcf, 'Fig4a_Bob_BER_multi_alpha.eps', '-depsc', '-r300');
print(gcf, 'Fig4a_Bob_BER_multi_alpha.pdf', '-dpdf', '-r300');

%% ================================================================
%  SECTION 3: FIG. 4(b) - 3D Surface Plot
%             Shows BER as a function of both delta and alpha for 
%             Bob (blue surface) and Eve (red surface).
%             Demonstrates the security gap across all parameter combinations.
% ================================================================
fprintf('\n=== Fig. 4(b): 3D BER Surface (Bob vs Eve) ===\n');

% Preallocate 3D BER arrays
BER_Bob_3D = zeros(numel(alpha_3D), numel(delta_dB_3D));
BER_Eve_3D = zeros(numel(alpha_3D), numel(delta_dB_3D));

for ia = 1:numel(alpha_3D)
    alpha = alpha_3D(ia);
    
    % Variances (Eq. (2) style)
    sigmaL2 = 2*P/(1+alpha);
    sigmaH2 = alpha*sigmaL2;
    
    % Generate bit stream for this alpha
    bits = randi([0 1], 1, nbits_3D);
    
    fprintf('  alpha = %.2f ...\n', alpha);
    
    for iD = 1:numel(delta_dB_3D)
        delta_dB = delta_dB_3D(iD);
        delta = 10^(delta_dB/10);
        sigmaW2 = sigmaL2 / delta;
        
        errB = 0; errE = 0;
        
        for k = 1:nbits_3D
            b = bits(k);
            
            % Channels (Rayleigh fading)
            h  = (randn + 1j*randn)/sqrt(2);
            he = (randn + 1j*randn)/sqrt(2);
            
            rho  = abs(h);
            rhoe = abs(he);
            
            % NoiseMod symbol generation
            if b == 0
                x = sqrt(sigmaL2) * randn(1, N);
            else
                x = sqrt(sigmaH2) * randn(1, N);
            end
            
            % CAS-NoiseMod pre-scaling: s[n] = rho * x[n]
            s_CAS = rho * x;
            
            % Independent noises for Bob and Eve
            wB = sqrt(sigmaW2/2) * (randn(1, N) + 1j*randn(1, N));
            wE = sqrt(sigmaW2/2) * (randn(1, N) + 1j*randn(1, N));
            
            % Received signals
            yB = h  * s_CAS + wB;
            yE = he * s_CAS + wE;
            
            % Decision statistics
            ShatB = mean(abs(yB).^2);
            ShatE = mean(abs(yE).^2);
            
            % Bob's matched threshold (knows rho)
            C0B = (rho^4) * sigmaL2 + sigmaW2;
            C1B = (rho^4) * sigmaH2 + sigmaW2;
            gammaB = (C1B*C0B)/(C1B-C0B) * log(C1B/C0B);
            bhatB = (ShatB > gammaB);
            
            % Eve's mismatched threshold (uses her own channel)
            C0E = (rhoe^4) * sigmaL2 + sigmaW2;
            C1E = (rhoe^4) * sigmaH2 + sigmaW2;
            gammaE = (C1E*C0E)/(C1E-C0E) * log(C1E/C0E);
            bhatE = (ShatE > gammaE);
            
            errB = errB + xor(bhatB, b);
            errE = errE + xor(bhatE, b);
        end
        
        BER_Bob_3D(ia, iD) = errB / nbits_3D;
        BER_Eve_3D(ia, iD) = errE / nbits_3D;
    end
end

%% Plot Fig. 4(b) - 3D Surface Plot showing Bob and Eve BER
[Delta_mesh, Alpha_mesh] = meshgrid(delta_dB_3D, alpha_3D);

figure('Color','white','Position',[100 100 900 700]);

% Plot Bob's BER surface (blue)
surf(Delta_mesh, Alpha_mesh, log10(BER_Bob_3D), ...
    'FaceColor', 'b', 'FaceAlpha', 0.7, 'EdgeColor', 'b', 'EdgeAlpha', 0.3);
hold on;

% Plot Eve's BER surface (red)
surf(Delta_mesh, Alpha_mesh, log10(BER_Eve_3D), ...
    'FaceColor', 'r', 'FaceAlpha', 0.7, 'EdgeColor', 'r', 'EdgeAlpha', 0.3);

% Axis labels and formatting
xlabel('\delta (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('\alpha', 'FontSize', 14, 'FontName', 'Times New Roman');
zlabel('BER', 'FontSize', 14, 'FontName', 'Times New Roman');

% Set Z-axis tick labels to show actual BER values
zticks = [-3 -2 -1 0];
set(gca, 'ZTick', zticks);
set(gca, 'ZTickLabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'});
zlim([-3 0]);

% View angle for better visualization
view(45, 25);

% Grid and box
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
set(gca, 'Box', 'on');

% Legend
legend({'Bob', 'Eve'}, 'Location', 'northeast', 'FontSize', 12);
title('(b)', 'FontSize', 14, 'FontName', 'Times New Roman');

% Colorbar (optional, shows BER gradient)
colormap([linspace(0,0,64)' linspace(0,0,64)' linspace(0.5,1,64)'; ...
          linspace(0.5,1,64)' linspace(0,0,64)' linspace(0,0,64)']);

print(gcf, 'Fig4b_3D_BER_Surface.eps', '-depsc', '-r300');
print(gcf, 'Fig4b_3D_BER_Surface.pdf', '-dpdf', '-r300');
print(gcf, 'Fig4b_3D_BER_Surface.png', '-dpng', '-r300');

%% ================================================================
%  SECTION 4: FIG. 4(c) - 2D Comparison Plot
%             Bob and Eve BER for CAS-NoiseMod vs Conventional NoiseMod
%             at fixed alpha. Demonstrates that the security gap vanishes
%             without channel-adaptive pre-scaling.
% ================================================================
fprintf('\n=== Fig. 4(c): 2D Comparison (CAS vs Baseline) at alpha = %.2f ===\n', alpha_F5);

% Preallocate BER arrays
BER_Bob_CAS  = zeros(1, numel(delta_F5));
BER_Eve_CAS  = zeros(1, numel(delta_F5));
BER_Bob_BASE = zeros(1, numel(delta_F5));
BER_Eve_BASE = zeros(1, numel(delta_F5));

% Variances for fixed alpha
sigmaL2 = 2*P/(1+alpha_F5);
sigmaH2 = alpha_F5*sigmaL2;

% Generate bit stream
bits = randi([0 1], 1, nbits_F5);

for iD = 1:numel(delta_F5)
    delta = delta_F5(iD);
    sigmaW2 = sigmaL2 / delta;

    errB_CAS = 0; errE_CAS = 0;
    errB_BAS = 0; errE_BAS = 0;

    for k = 1:nbits_F5
        b = bits(k);

        % Channels (Rayleigh fading)
        h  = (randn + 1j*randn)/sqrt(2);
        he = (randn + 1j*randn)/sqrt(2);

        rho  = abs(h);
        rhoe = abs(he);

        % Symbol generation: x[n] real Gaussian
        if b == 0
            x = sqrt(sigmaL2) * randn(1, N);
        else
            x = sqrt(sigmaH2) * randn(1, N);
        end

        % Independent noises
        wB = sqrt(sigmaW2/2) * (randn(1, N) + 1j*randn(1, N));
        wE = sqrt(sigmaW2/2) * (randn(1, N) + 1j*randn(1, N));

        % --------------------------
        % (A) CAS-NoiseMod (with channel-adaptive pre-scaling)
        % --------------------------
        s_CAS = rho * x;               % pre-scaling by channel magnitude
        yB_CAS = h  * s_CAS + wB;
        yE_CAS = he * s_CAS + wE;

        ShatB = mean(abs(yB_CAS).^2);
        ShatE = mean(abs(yE_CAS).^2);

        % Bob matched threshold under CAS: Cb = rho^4*sigma_b^2 + sigmaW2
        C0B = (rho^4) * sigmaL2 + sigmaW2;
        C1B = (rho^4) * sigmaH2 + sigmaW2;
        gammaB = (C1B*C0B)/(C1B-C0B) * log(C1B/C0B);
        bhatB = (ShatB > gammaB);

        % Eve mismatched threshold under CAS (uses her own channel)
        C0E_assumed = (rhoe^4) * sigmaL2 + sigmaW2;
        C1E_assumed = (rhoe^4) * sigmaH2 + sigmaW2;
        gammaE_mis  = (C1E_assumed*C0E_assumed)/(C1E_assumed-C0E_assumed) * log(C1E_assumed/C0E_assumed);
        bhatE_CAS = (ShatE > gammaE_mis);

        errB_CAS = errB_CAS + xor(bhatB, b);
        errE_CAS = errE_CAS + xor(bhatE_CAS, b);

        % --------------------------
        % (B) Conventional NoiseMod baseline (no pre-scaling)
        % --------------------------
        if INCLUDE_BASELINE_NOISEMOD
            s_BASE  = x;              % no pre-scaling
            yB_BASE = h  * s_BASE + wB;
            yE_BASE = he * s_BASE + wE;

            ShatB0 = mean(abs(yB_BASE).^2);
            ShatE0 = mean(abs(yE_BASE).^2);

            % Bob threshold baseline: Cb = |h|^2*sigma_b^2 + sigmaW2
            C0B0 = (abs(h)^2)  * sigmaL2 + sigmaW2;
            C1B0 = (abs(h)^2)  * sigmaH2 + sigmaW2;
            gammaB0 = (C1B0*C0B0)/(C1B0-C0B0) * log(C1B0/C0B0);
            bhatB0 = (ShatB0 > gammaB0);

            % Eve threshold baseline: matched to her own channel (she knows he)
            C0E0 = (abs(he)^2) * sigmaL2 + sigmaW2;
            C1E0 = (abs(he)^2) * sigmaH2 + sigmaW2;
            gammaE0 = (C1E0*C0E0)/(C1E0-C0E0) * log(C1E0/C0E0);
            bhatE0 = (ShatE0 > gammaE0);

            errB_BAS = errB_BAS + xor(bhatB0, b);
            errE_BAS = errE_BAS + xor(bhatE0, b);
        end
    end

    BER_Bob_CAS(iD)  = errB_CAS / nbits_F5;
    BER_Eve_CAS(iD)  = errE_CAS / nbits_F5;
    BER_Bob_BASE(iD) = errB_BAS / nbits_F5;
    BER_Eve_BASE(iD) = errE_BAS / nbits_F5;

    fprintf('  delta = %d dB done\n', delta_dB_F5(iD));
end

%% Plot Fig. 4(c) - 2D Comparison: CAS-NoiseMod vs Conventional NoiseMod
figure('Color','white','Position',[100 100 850 600]);
hold on; grid on;
set(gca, 'YScale', 'log');

% Plot CAS-NoiseMod curves (proposed scheme)
h1 = semilogy(delta_dB_F5, BER_Bob_CAS, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerIndices', 1:2:length(delta_dB_F5), ...
    'DisplayName', 'Bob (CAS-NoiseMod)');
h2 = semilogy(delta_dB_F5, BER_Eve_CAS, '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerIndices', 1:2:length(delta_dB_F5), ...
    'DisplayName', 'Eve (CAS-NoiseMod)');

% Plot Baseline NoiseMod curves (conventional scheme without pre-scaling)
h4 = semilogy(delta_dB_F5, BER_Eve_BASE, '-d', 'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerIndices', 2:2:length(delta_dB_F5), ...
    'DisplayName', 'Eve (NoiseMod)');
h3 = semilogy(delta_dB_F5, BER_Bob_BASE, '-^', 'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerIndices', 2:2:length(delta_dB_F5), ...
    'DisplayName', 'Bob (NoiseMod)');

% Axis formatting
xlabel('\delta (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('BER', 'FontSize', 14, 'FontName', 'Times New Roman');
xlim([min(delta_dB_F5) max(delta_dB_F5)]);
ylim([1e-5 1]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');

% Black frame around the plot
set(gca, 'Box', 'on', 'LineWidth', 1.5);

% Legend in northeast corner
lg = legend([h1, h2, h3, h4], 'Location', 'northeast');
set(lg, 'FontSize', 11, 'Box', 'on');
title('(c)', 'FontSize', 14, 'FontName', 'Times New Roman');

% Save figures
print(gcf, 'Fig4c_2D_CAS_vs_Baseline.eps', '-depsc', '-r300');
print(gcf, 'Fig4c_2D_CAS_vs_Baseline.pdf', '-dpdf', '-r300');
print(gcf, 'Fig4c_2D_CAS_vs_Baseline.png', '-dpng', '-r300');

fprintf('\n=== All figures generated successfully ===\n');

%% ================================================================
%  FUNCTION: Bob theoretical BEP for CAS-NoiseMod (moment-matched Gamma)
%  Implements the 1D averaging over r = |h|^2 with exp(-r) pdf.
%  Matches the paper's Sec. III flow (Eq. (12)-(21) structure).
% ================================================================
function BEP = BEP_Bob_CAS_Theory(delta_dB_vec, N, alpha, P)

    delta_vec = 10.^(delta_dB_vec/10);

    % Variances (paper Eq. (2) form)
    sigmaL2 = 2*P/(1+alpha);
    sigmaH2 = alpha*sigmaL2;

    BEP = zeros(size(delta_vec));

    for iD = 1:numel(delta_vec)
        delta  = delta_vec(iD);
        sigmaW2 = sigmaL2 / delta;

        % Use r = |h|^2 ~ Exp(1). Under CAS, conditional means:
        % Cb(r) = r^2 * sigma_b^2 + sigmaW2  (dual-channel effect)
        C0 = @(r) (r.^2) * sigmaL2 + sigmaW2;
        C1 = @(r) (r.^2) * sigmaH2 + sigmaW2;

        % Impropriety factor (paper Eq. (13) form): zeta_b = ((Cb - sigmaW2)^2)/(Cb^2)
        zeta0 = @(r) ((C0(r) - sigmaW2).^2) ./ (C0(r).^2);
        zeta1 = @(r) ((C1(r) - sigmaW2).^2) ./ (C1(r).^2);

        % Moment-matched Gamma parameters (paper Eq. (14) form)
        k0 = @(r) N ./ (1 + zeta0(r));
        k1 = @(r) N ./ (1 + zeta1(r));

        th0 = @(r) C0(r) .* (1 + zeta0(r)) ./ N;
        th1 = @(r) C1(r) .* (1 + zeta1(r)) ./ N;

        % Threshold: use equal-shape approximation style (paper Eq. (16))
        gamma = @(r) (C1(r).*C0(r)) ./ (C1(r)-C0(r)) .* log(C1(r)./C0(r));

        % Error probabilities under Gamma approximation:
        % P10 = 1 - F_Gamma(gamma; k0, th0)
        % P01 = F_Gamma(gamma; k1, th1)
        integrand = @(r) 0.5 .* ( ...
            (1 - gamcdf(gamma(r), k0(r), th0(r))) + ...
            (    gamcdf(gamma(r), k1(r), th1(r))) ) .* exp(-r);

        BEP(iD) = integral(integrand, 0, Inf, ...
            'RelTol', 1e-7, 'AbsTol', 1e-12, 'ArrayValued', true);
    end
end

