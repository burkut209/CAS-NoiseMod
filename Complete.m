%% ================================================================
%  CAS-NoiseMod (WCL2025-2845) - UNIFIED MATLAB SCRIPT
%  Author: Serdar Prencuva (clean unified version)
%
%  Purpose (matches paper):
%   (A) CAS-NoiseMod simulation: Bob BER, Eve BER (mismatched threshold)
%   (B) Bob theoretical BEP via moment-matched Gamma averaging
%   (C) Reviewer-required baseline: Conventional NoiseMod (no pre-scaling)
%       - Bob BER baseline
%       - Eve BER baseline
%   (D) Figures kept in the same spirit as the paper:
%       - Fig. 4: Bob BER vs delta for multiple alpha (simulation + theory)
%       - Fig. 5: 3D surfaces (delta, alpha) for Bob vs Eve
%         with baseline overlaid (still ONE figure, no extra figures)
%
%  Notes:
%   - CAS-NoiseMod model: s[n] = rho * x[n], rho = |h| (dual-channel effect),
%     y[n] = h*s[n] + w[n]. See System Model Sec. II.  (Eq. (3)-(6)) 
%   - Conventional NoiseMod baseline: s[n] = x[n], y[n] = h*x[n] + w[n]
%   - Eve in CAS-NoiseMod: y_e[n] = h_e*s[n] + w_e[n],
%     but Eve does not know rho = |h| and thus uses a mismatched threshold.
%   - Random seed fixed for repeatability (reviewer request).
% ================================================================

clear; close all; clc;

%% ================================================================
%  SECTION 1: GLOBAL PARAMETERS (paper-consistent knobs)
% ================================================================
rng(3, 'twister');                  % Fixed seed for repeatability

N        = 120;                     % Samples per bit (paper uses N=120)
P        = 1;                       % Average transmit power
nbits_F4 = 1e5;                     % Bits for Fig. 4 style curves (heavy but consistent)
nbits_F5 = 2e4;                     % Bits per (alpha,delta) point for Fig. 5 surface (increase if needed)

delta_dB_F4 = -10:1:35;             % delta(dB) sweep for Fig. 4
delta_F4    = 10.^(delta_dB_F4/10);

alpha_list_F4 = [1, 1.5, 3, 10];    % Fig. 4 alpha set (as used in your current plotting)

% Fig. 5 surface resolution (delta vs alpha)
delta_dB_F5 = -10:2:30;
delta_F5    = 10.^(delta_dB_F5/10);

alpha_grid_F5 = 1:0.5:10;           % alpha range (paper says 1 to 10)

% Toggle baseline overlay (reviewer required)
INCLUDE_BASELINE_NOISEMOD = true;

%% ================================================================
%  SECTION 2: FIG. 4 (Bob BER vs delta for multiple alpha)
%            Simulation + Theoretical (CAS-NoiseMod / Bob only)
%            (This section reproduces the role of Fig. 4 in the paper.)
% ================================================================
BER_Bob_sim_F4  = zeros(numel(alpha_list_F4), numel(delta_F4));
BER_Bob_the_F4  = zeros(numel(alpha_list_F4), numel(delta_F4));

fprintf('=== Fig. 4: Bob BER (CAS-NoiseMod) for multiple alpha ===\n');

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

            % Threshold (equal-shape / variance-threshold form; Eq. (8)/(17) spirit)
            % Under CAS, Bob’s conditional variances: Cb = (rho^4)*sigma_b^2 + sigmaW2
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

%% Plot Fig. 4 (keep ONE figure)
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

xlabel('\delta (dB)');
ylabel('BER');
ylim([1e-3 1]);
xlim([min(delta_dB_F4) max(delta_dB_F4)]);
legend('Location','northeast');

print(gcf, 'Fig4_Bob_BER_multi_alpha.eps', '-depsc', '-r300');
print(gcf, 'Fig4_Bob_BER_multi_alpha.pdf', '-dpdf', '-r300');

%% ================================================================
%  SECTION 3: FIG. 5 (3D surfaces: Bob vs Eve, CAS-NoiseMod)
%            + Reviewer baseline overlay: conventional NoiseMod
%            (still ONE figure)
% ================================================================
fprintf('\n=== Fig. 5: 3D surfaces (CAS-NoiseMod) + baseline overlay ===\n');

BER_Bob_CAS_F5 = zeros(numel(alpha_grid_F5), numel(delta_F5));
BER_Eve_CAS_F5 = zeros(numel(alpha_grid_F5), numel(delta_F5));

BER_Bob_BASE_F5 = zeros(numel(alpha_grid_F5), numel(delta_F5));
BER_Eve_BASE_F5 = zeros(numel(alpha_grid_F5), numel(delta_F5));

for ia = 1:numel(alpha_grid_F5)
    alpha = alpha_grid_F5(ia);

    sigmaL2 = 2*P/(1+alpha);
    sigmaH2 = alpha*sigmaL2;

    bits = randi([0 1], 1, nbits_F5);

    for iD = 1:numel(delta_F5)
        delta = delta_F5(iD);
        sigmaW2 = sigmaL2 / delta;

        errB_CAS = 0; errE_CAS = 0;
        errB_BAS = 0; errE_BAS = 0;

        for k = 1:nbits_F5
            b = bits(k);

            % Channels
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

            % Independent noises (more correct than reusing same noise)
            wB = sqrt(sigmaW2/2) * (randn(1, N) + 1j*randn(1, N));
            wE = sqrt(sigmaW2/2) * (randn(1, N) + 1j*randn(1, N));

            % --------------------------
            % (A) CAS-NoiseMod
            % --------------------------
            s_CAS = rho * x;               % pre-scaling
            yB_CAS = h  * s_CAS + wB;
            yE_CAS = he * s_CAS + wE;

            ShatB = mean(abs(yB_CAS).^2);
            ShatE = mean(abs(yE_CAS).^2);

            % Bob matched threshold under CAS: Cb = rho^4*sigma_b^2 + sigmaW2
            C0B = (rho^4) * sigmaL2 + sigmaW2;
            C1B = (rho^4) * sigmaH2 + sigmaW2;
            gammaB = (C1B*C0B)/(C1B-C0B) * log(C1B/C0B);
            bhatB = (ShatB > gammaB);

            % Eve mismatched threshold under CAS:
            % Actual Eve variance depends on (|he|^2)*(|h|^2)*sigma_b^2 + sigmaW2,
            % but Eve does not know |h|. A standard mismatch model is:
            %   Eve assumes the pre-scaling used rho_assumed = |he|
            % which yields assumed Cb,e = |he|^4*sigma_b^2 + sigmaW2.
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

        BER_Bob_CAS_F5(ia, iD) = errB_CAS / nbits_F5;
        BER_Eve_CAS_F5(ia, iD) = errE_CAS / nbits_F5;

        if INCLUDE_BASELINE_NOISEMOD
            BER_Bob_BASE_F5(ia, iD) = errB_BAS / nbits_F5;
            BER_Eve_BASE_F5(ia, iD) = errE_BAS / nbits_F5;
        end
    end
end

%% --- Figure 5 Plotting Block (REPLACE your existing Fig. 5 plotting) ---

figure('Color','white','Position',[150 120 950 700], 'Renderer','opengl');
hold on;

% Meshgrid for Fig.5 vectors
[X, Y] = meshgrid(delta_dB_F5, alpha_grid_F5);

% Clamp to zmin to avoid ugly spikes below axis floor
zmin = 1e-3;
Z_Bob_CAS  = max(BER_Bob_CAS_F5 , zmin);
Z_Eve_CAS  = max(BER_Eve_CAS_F5 , zmin);
Z_Bob_BASE = max(BER_Bob_BASE_F5, zmin);
Z_Eve_BASE = max(BER_Eve_BASE_F5, zmin);

% =====================================================================
% COLOR DEFINITIONS (clean, distinct)
% =====================================================================
% CAS colors
bobCAS_Blue = [0 0.4470 0.7410];      % MATLAB blue
eveCAS_Red  = [0.8 0.1 0.1];          % Deep red (much redder!)

% Baseline colors (clearly distinct)
bobBASE_Cyan   = [0.3010 0.7450 0.9330]; % Light cyan/sky blue
eveBASE_Purple = [0.4940 0.1840 0.5560]; % MATLAB purple

% =====================================================================
% CAS-NoiseMod Surfaces
% =====================================================================
hBobCAS = surf(X, Y, Z_Bob_CAS, ...
    'FaceColor', bobCAS_Blue, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hEveCAS = surf(X, Y, Z_Eve_CAS, ...
    'FaceColor', eveCAS_Red, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% =====================================================================
% NoiseMod Baseline Surfaces
% =====================================================================
hBobBASE = surf(X, Y, Z_Bob_BASE, ...
    'FaceColor', bobBASE_Cyan, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hEveBASE = surf(X, Y, Z_Eve_BASE, ...
    'FaceColor', eveBASE_Purple, 'EdgeColor', 'none', 'FaceAlpha', 0.4);

% =====================================================================
% Axes formatting
% =====================================================================
set(gca, 'ZScale', 'log');
zlim([1e-3 1]);
set(gca, 'ZTick', [1e-3 1e-2 1e-1 1e0]);

xlabel('\delta (dB)', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('\alpha',      'FontSize', 16, 'FontName', 'Times New Roman');
zlabel('BER',         'FontSize', 16, 'FontName', 'Times New Roman');

xlim([min(delta_dB_F5) max(delta_dB_F5)]);
ylim([min(alpha_grid_F5) max(alpha_grid_F5)]);

view(52, 22);

grid on;
set(gca, 'GridAlpha', 0.12, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out');
box on;

% Disable lighting to keep colors clean/solid
lighting none;

lg = legend([hBobCAS, hEveCAS, hBobBASE, hEveBASE], ...
    {'Bob (CAS)', 'Eve (CAS, mismatched)', 'Bob (NoiseMod baseline)', 'Eve (NoiseMod baseline)'}, ...
    'Location','northwest');
set(lg, 'FontSize', 12, 'Box', 'on');

print(gcf, '3d_ber_performance.eps', '-depsc', '-r300');
print(gcf, '3d_ber_performance.pdf', '-dpdf',  '-r300');
print(gcf, '3d_ber_performance.png', '-dpng',  '-r300');

%% ================================================================
%  FUNCTION: Bob theoretical BEP for CAS-NoiseMod (moment-matched Gamma)
%  Implements the 1D averaging over r = |h|^2 with exp(-r) pdf.
%  Matches the paper’s Sec. III flow (Eq. (12)-(21) structure).
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

        % Moment-matched Gamma parameters (paper Eq. (15) form)
        k0 = @(r) N ./ (1 + zeta0(r));
        k1 = @(r) N ./ (1 + zeta1(r));

        th0 = @(r) C0(r) .* (1 + zeta0(r)) ./ N;
        th1 = @(r) C1(r) .* (1 + zeta1(r)) ./ N;

        % Threshold: use equal-shape approximation style (paper Eq. (17))
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
