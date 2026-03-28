% Unified optical-guided phase-only calibration
% Supports: MFA and RMA
%
% Dataset convention:
%   1 -> YANIK
%   2 -> 3DRIED case 1
%   3 -> 3DRIED case 2
%   4 -> OURS case 1
%   5 -> OURS case 2

%% ============================================================
% USER SELECTION
% ============================================================
clc; clear; close all;

sar_algo = 'MFA';     % 'MFA' or 'RMA'
dataset  = 3;         % 1..5

% Directories
optical_ref_dir = fullfile(pwd, 'optical_refs');
raw_dataDir     = fullfile(pwd, 'sar_raw_data');

% Constants
c0 = physconst('lightspeed');

% ============================================================
% DATASET-SPECIFIC PREPROCESSING
% ============================================================

if dataset == 1
    % ==============================================================
    % DATASET 1 YANIK
    % Note: TrueImage size 76x306 ---> ref_optical_image 76x306
    % ==============================================================
    obj_name = "yanik";
    dataName = "rawData3D_simple2D.mat";

    bbox = [-70 80 -80 70];
    z0   = 280;                        % mm

    nFFTtime  = 1024;
    nFFTspace = 1024;

    F0 = 77e9;
    FS = 9121e3;
    K0 = 63.343e12;
    tI = 4.5225e-10;

    % Optical reference for dataset 1
    optical_ref_file    = 'yanik_ref_optical_image.mat';
    optical_ref_varname = 'ref_optical_image';
    invert_optical_reference = true;

    rawStruct = load(fullfile(raw_dataDir, dataName));
    if isfield(rawStruct, 'rawData3D_simple2D')
        sarRawData = rawStruct.rawData3D_simple2D;
    else
        fn = fieldnames(rawStruct);
        sarRawData = rawStruct.(fn{1});
    end

    [Nsamp, M, N] = size(sarRawData);

    dx = 200/(N-1);
    dy = 198/(M-1);

elseif dataset == 2 || dataset == 3
    % ==============================================================
    % DATASET 2 or 3 (3DRIED)
    % Note: TrueImage size 201x401 ---> ref_optical_image 201x401
    % ==============================================================
    obj_name = "3dried";

    if dataset == 2
        dataName = "rawdata_gun_3dried.mat";
        optical_ref_file = 'ref_optical_image_gun_3dried.mat';
        bbox = [-150 250 -320 80];
    else
        dataName = "rawdata_knife_3dried.mat";
        optical_ref_file = 'ref_optical_image_knife_3dried.mat';
        bbox = [-160 250 -210 260];
    end

    optical_ref_varname = 'ref_optical_image';
    invert_optical_reference = true;

    %bbox = [-150 250 -320 80];

    nFFTspace = 512;
    F0 = (77 + 1.8)*1e9;
    FS = 5e6;
    K0 = 70.295e12;
    tI = 6.2516e-10;

    dx = 1;
    dy = 2;

    rawStruct = load(fullfile(raw_dataDir, dataName));

    if isfield(rawStruct, 'Echo')
        Echo = rawStruct.Echo;
    elseif isfield(rawStruct, 'adcRawData')
        adcRawData = rawStruct.adcRawData;

        m = length(adcRawData.data);
        num_TX = 1;
        num_RX = 4;
        num_sample = 256;

        Raw_echo = zeros(m*num_TX*num_RX, num_sample);
        for ii = 1:num_TX*m
            Raw_echo((ii-1)*4+1:ii*4, :) = squeeze(adcRawData.data{ii});
        end

        Num = num_TX*num_RX;
        Nx  = 407;
        Nz  = 200;

        Sr0  = Raw_echo(1:Num:end, :);
        Echo = zeros(Nx*Nz, num_sample);

        err = 43;
        for ii = 1:Nz
            kk = floor(ii/Nz*err);
            Echo((ii-1)*Nx+1:ii*Nx, :) = Sr0((ii-1)*Nx+1+kk:ii*Nx+kk, :);
        end
        Echo = reshape(Echo, [Nx, Nz, num_sample]);
    else
        error('Dataset 2/3 file must contain Echo or adcRawData.');
    end

    [Nx, Nz, num_sample] = size(Echo);

    sarRawData = Echo;
    Nsamp = num_sample;
    nFFTtime = num_sample;
    M = Nz;
    N = Nx;

    ID_select = 17;
    z0 = (c0/2) * ((ID_select/(K0*(1/FS)*Nsamp)) - tI) * 1000;   % mm

elseif dataset == 4 || dataset == 5
    % ==============================================================
    % DATASET 4 or 5 OURS
    % ==============================================================
    obj_name = "ours";

    if dataset == 4
        dataName = "rawdata_knife_ours.mat";
        optical_ref_file = 'ref_optical_image_knife_ours.mat';
        dx = 1; dy = 1;
        z0   = 185;  
        bbox = [-200 200 -200 200]; 
    else
        dataName = "rawdata_plier_ours.mat";
        optical_ref_file = 'ref_optical_image_plier_ours.mat';
        dx = 1; dy = 2;
        z0   = 210; 
        bbox = [-200 200 -200 200];  
    end

    optical_ref_varname = 'ref_optical_image';
    invert_optical_reference = true;

    %bbox = [-70 80 -80 70];
    %z0   = 280;

    nFFTtime  = 1024;
    nFFTspace = 1024;

    F0 = 77e9;
    FS = 5000e3;
    K0 = 70.295e12;
    tI = 4.5225e-10;

    rawStruct = load(fullfile(raw_dataDir, dataName));
    if isfield(rawStruct, 'adcDataCube')
        sarRawData = rawStruct.adcDataCube;
    else
        fn = fieldnames(rawStruct);
        sarRawData = rawStruct.(fn{1});
    end
    %for ii = 2:2:size(sarRawData, 1)
    %    sarRawData(ii, :) = fliplr(sarRawData(ii, :));
    %end
    [Nsamp, M, N] = size(sarRawData);


else
    error('Unknown dataset selection. Choose dataset = 1, 2, 3, 4, or 5.');
end

fprintf('Selected dataset: %d\n', dataset);
fprintf('Selected family : %s\n', upper(obj_name));
fprintf('Selected algorithm: %s\n', upper(sar_algo));
fprintf('dx = %.3f mm, dy = %.3f mm, z0 = %.3f mm, FS = %.3f MHz\n', ...
    dx, dy, z0, FS/1e6);

% ============================================================
% LOAD OPTICAL REFERENCE
% ============================================================

opticalStruct = load(fullfile(optical_ref_dir, optical_ref_file));
if ~isfield(opticalStruct, optical_ref_varname)
    error('Variable "%s" not found in file "%s".', ...
        optical_ref_varname, optical_ref_file);
end

ref_optical_image = double(opticalStruct.(optical_ref_varname));

if invert_optical_reference
    ref_optical_image = 1 - ref_optical_image;
end

if any(dataset == [1 2 3 4 5])
    ref_optical_image = fliplr(ref_optical_image);
end

% ============================================================
% PARAMETER STRUCT
% ============================================================

params = struct;
params.sar_algo   = upper(sar_algo);
params.bbox       = bbox;
params.dx         = dx;
params.dy         = dy;
params.z0         = z0;            % mm for MFA; overwritten in meters for RMA
params.F0         = F0;
params.FS         = FS;
params.K0         = K0;
params.tI         = tI;
params.nFFTspace  = nFFTspace;
params.Nsamp      = Nsamp;
params.M          = M;
params.N          = N;

if exist('nFFTtime', 'var')
    params.nFFTtime = nFFTtime;
end

%%============================================================
% BASELINE RECONSTRUCTION SETUP
% ============================================================

switch upper(sar_algo)

    case 'MFA'
        %if dataset == 2 || dataset == 3
        %    error('MFA branch currently expects rawData3D-style input. Dataset 2/3 is currently set up for 3DRIED/RMA-style Echo data.');
        %end
        if dataset == 2 || dataset == 3
            rawDataFFT = fft(sarRawData, nFFTtime, 3);
    
            k0_range_bin = 17;
            sarData = squeeze(rawDataFFT(:, :, k0_range_bin)).';
    
            for ii = 2:2:size(sarData,1)
                sarData(ii,:) = fliplr(sarData(ii,:));
            end
    
            params.k0_range_bin = k0_range_bin;
         
        else
            k0_range_bin = round(K0 / FS * (2 * z0 * 1e-3 / c0 + tI) * nFFTtime);
 
            rawDataFFT = fft(sarRawData, nFFTtime, 1);
            sarData    = squeeze(rawDataFFT(k0_range_bin + 1, :, :));
            
            if dataset == 4 || dataset == 5
                for ii = 2:2:size(sarData,1)
                    sarData(ii,:) = fliplr(sarData(ii,:));
                end
            end
            %sarData = sarData(:, end:-1:1);
            params.k0_range_bin = k0_range_bin;
        end

        %k0_range_bin = round(K0 / FS * (2 * z0 * 1e-3 / c0 + tI) * nFFTtime);
        %rawDataFFT = fft(sarRawData, nFFTtime, 1);
        %sarData    = squeeze(rawDataFFT(k0_range_bin + 1, :, :));

        %sarData = sarData(:, end:-1:1);
        %params.k0_range_bin = k0_range_bin;

        H_ref = refMF(params);
        params.operator_ref   = H_ref;
        params.operator_abs   = abs(H_ref);
        params.operator_phase = angle(H_ref);

        [xAxis, yAxis, baseline_img_dl, ~] = dlMFA_custom(sarData, H_ref, params);
        baseline_img = extractdata(baseline_img_dl);

    case 'RMA'
        if dataset == 2 || dataset == 3
            Echo = sarRawData;                      % Nx x Nz x Nsamp
            num_sample = size(Echo, 3);
            rawDataFFT = fft(Echo, num_sample, 3);

            k0_range_bin = 17;
            sarData = squeeze(rawDataFFT(:, :, k0_range_bin)).';   % Nz x Nx

            for ii = 2:2:size(sarData,1)
                sarData(ii,:) = fliplr(sarData(ii,:));
            end

            z0_eff = (c0/2) * (((k0_range_bin - 1) / (K0*(1/FS)*num_sample)) - tI);  % meters

            params.k0_range_bin = k0_range_bin;
            params.z0 = z0_eff;
            params.nFFTtime = num_sample;

        else
            Echo = permute(sarRawData, [3, 2, 1]);   % Nx x Ny x Nsamp
            [Nx, Nz, ~] = size(Echo);
            num_sample = size(Echo, 3);
            rawDataFFT = fft(Echo, num_sample, 3);

            E = squeeze(sum(sum(abs(rawDataFFT).^2, 1), 2));
            [~, k0_range_bin] = max(E);

            sarData = squeeze(rawDataFFT(:, :, k0_range_bin)).';
            %sarData = sarData(:, end:-1:1);
            if dataset == 4 || dataset == 5
                for ii = 2:2:Nz
                    sarData(ii, :) = fliplr(sarData(ii, :));
                end
            end

            z0_eff = (c0/2) * (((k0_range_bin - 1) / (K0*(1/FS)*num_sample)) - tI);  % meters

            params.k0_range_bin = k0_range_bin;
            params.z0 = z0_eff;
            params.nFFTtime = num_sample;
        end

        K_ref = buildRMAKernel(params);
        params.operator_ref   = K_ref;
        params.operator_abs   = abs(K_ref);
        params.operator_phase = angle(K_ref);

        [xAxis, yAxis, baseline_img_dl, ~] = dlRMA_custom(sarData, K_ref, params);
        baseline_img = extractdata(baseline_img_dl);

    otherwise
        error('sar_algo must be MFA or RMA.');
end

% ============================================================
% MATCH OPTICAL SIZE + GLOBAL NORMALIZATION
% ============================================================

if ~isequal(size(ref_optical_image), size(baseline_img))
    fprintf('Resizing optical reference %s -> %s\n', ...
        mat2str(size(ref_optical_image)), mat2str(size(baseline_img)));

    if exist('imresize','file') == 2
        ref_optical_image = imresize(ref_optical_image, size(baseline_img), 'bilinear');
    else
        [X, Y]  = meshgrid(linspace(0,1,size(ref_optical_image,2)), ...
                           linspace(0,1,size(ref_optical_image,1)));
        [Xq,Yq] = meshgrid(linspace(0,1,size(baseline_img,2)), ...
                           linspace(0,1,size(baseline_img,1)));
        ref_optical_image = interp2(X, Y, ref_optical_image, Xq, Yq, 'linear', 0);
    end
end

global_scale = max(abs(baseline_img(:))) + 1e-12;
params.global_scale = global_scale;

baseline_img      = baseline_img / global_scale;
ref_optical_image = ref_optical_image / (max(abs(ref_optical_image(:))) + 1e-12);

params.ref_optical_image = dlarray(ref_optical_image, "SS");

% ============================================================
% SHOW BASELINE + OPTICAL
% ============================================================

figure;
subplot(1,2,1);
imagesc(xAxis, yAxis, baseline_img);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca,'jet');
colorbar;
title(sprintf('Baseline %s image', upper(sar_algo)));

subplot(1,2,2);
imagesc(xAxis, yAxis, ref_optical_image);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca,'gray');
colorbar;
title('Optical reference');

%% ============================================================
% PHASE-ONLY OPTIMIZATION SETUP
% ============================================================

switch upper(sar_algo)
    case 'MFA'
        params.lambda_phi    = 1e-3;
        params.lambda_smooth = 1e-3;
        params.lambda_dev    = 1e-3;

        lr          = 1e3;
        maxIter_cap = 1000;

    case 'RMA'
        params.lambda_phi    = 1e-3;
        params.lambda_smooth = 1e-3;
        params.lambda_dev    = 1e-3;

        lr          = 7e3;
        maxIter_cap = 1000;

    otherwise
        error('sar_algo must be MFA or RMA.');
end

PhaseVar = dlarray(zeros(size(params.operator_ref), 'like', real(params.operator_ref)));

use_early_stop = true;
check_every   = 50;
W             = 30;
rel_drop_min  = 1e-3;
abs_drop_min  = 0;
K_weak        = 3;

iter = 0;
bestLoss = inf;
bestIter = 0;
bestPhaseVar = PhaseVar;
loss_hist = nan(maxIter_cap,1);
weakWinCount = 0;

fprintf(['\n=== Starting phase-only %s calibration | lr=%.3e | maxIter=%d | ' ...
         'lambda_phi=%.3e | lambda_smooth=%.3e | lambda_dev=%.3e ===\n'], ...
    upper(sar_algo), lr, maxIter_cap, ...
    params.lambda_phi, params.lambda_smooth, params.lambda_dev);

while true
    if iter >= maxIter_cap
        fprintf('Reached maxIter cap (%d).\n', maxIter_cap);
        break;
    end

    iter = iter + 1;

    [loss, grad, current_img] = dlfeval(@loss_and_grad_phaseOnly_MFA_RMA, sarData, PhaseVar, params);

    lossVal         = double(gather(extractdata(loss)));
    loss_hist(iter) = lossVal;

    PhaseVar = PhaseVar - lr * grad;

    if lossVal < bestLoss
        bestLoss     = lossVal;
        bestIter     = iter;
        bestPhaseVar = PhaseVar;
    end

    if mod(iter, 25) == 0 || iter == 1 || iter == maxIter_cap
        mse_img_vs_optical = double(gather(extractdata( ...
            mean((current_img - params.ref_optical_image).^2, 'all'))));

        grad_num = extractdata(grad);
        Gnow = max(abs(grad_num), [], 'all');

        fprintf(['Iter %04d | Loss=%.4e (best=%.4e @%d) | G=%.6e | ' ...
                 'MSE(image,optical)=%.4e\n'], ...
                iter, lossVal, bestLoss, bestIter, Gnow, mse_img_vs_optical);
    end

    if use_early_stop && iter > W && mod(iter, check_every) == 0
        L_old = loss_hist(iter - W);
        L_new = loss_hist(iter);

        rel_drop = (L_old - L_new) / max(abs(L_old), 1e-12);
        abs_drop = (L_old - L_new);

        if abs_drop_min > 0
            isWeak = (rel_drop < rel_drop_min) && (abs_drop < abs_drop_min);
        else
            isWeak = (rel_drop < rel_drop_min);
        end

        if isWeak
            weakWinCount = weakWinCount + 1;
        else
            weakWinCount = 0;
        end

        if weakWinCount >= K_weak
            fprintf(['Early stop: diminishing returns. Over last %d iters: ' ...
                     'rel_drop=%.3e, abs_drop=%.3e. Triggered %d/%d.\n'], ...
                     W, rel_drop, abs_drop, weakWinCount, K_weak);
            break;
        end
    end
end

loss_hist = loss_hist(1:iter);
PhaseVar = bestPhaseVar;

fprintf('-----------------------------\n');
fprintf('Done. Best loss %.4e at iter %d (ran %d iters).\n', bestLoss, bestIter, iter);
fprintf('-----------------------------\n');

%% ============================================================
% FINAL OPTIMIZED RECONSTRUCTION
% ============================================================

optimized_img_dl = reconstruct_phaseOnly_MFA_RMA(sarData, PhaseVar, params);
optimized_img = extractdata(optimized_img_dl);
optimized_img = optimized_img / params.global_scale;

figure;
imagesc(xAxis, yAxis, optimized_img);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca,'jet');
colorbar;
title(sprintf('Optimized %s image', upper(sar_algo)));

%% ============================================================
% NO-REFERENCE IMAGE QUALITY METRICS
% ============================================================

m_base = image_quality_metrics(baseline_img);
m_opt  = image_quality_metrics(optimized_img);

fprintf('\n================ Image quality summary (NO-REFERENCE) ================\n');
fprintf('Baseline %s:\n', upper(sar_algo));
fprintf('  Contrast : %.6f\n', m_base.contrast);
fprintf('  Entropy  : %.6f\n', m_base.entropy);
fprintf('  Peak     : %.6e\n', m_base.peak);

fprintf('\nOptimized %s:\n', upper(sar_algo));
fprintf('  Contrast : %.6f\n', m_opt.contrast);
fprintf('  Entropy  : %.6f\n', m_opt.entropy);
fprintf('  Peak     : %.6e\n', m_opt.peak);
fprintf('=======================================================\n');

%% ============================================================
% FINAL METRICS
% ============================================================

mse_opt_vs_optical  = mean((optimized_img(:) - ref_optical_image(:)).^2);
mse_opt_vs_baseline = mean((optimized_img(:) - baseline_img(:)).^2);
mse_base_vs_optical = mean((baseline_img(:)  - ref_optical_image(:)).^2);

num1 = sum(optimized_img(:) .* ref_optical_image(:));
den1 = sqrt(sum(optimized_img(:).^2) * sum(ref_optical_image(:).^2)) + 1e-12;
ncc_opt_vs_optical = num1 / den1;

num2 = sum(optimized_img(:) .* baseline_img(:));
den2 = sqrt(sum(optimized_img(:).^2) * sum(baseline_img(:).^2)) + 1e-12;
ncc_opt_vs_baseline = num2 / den2;

num3 = sum(baseline_img(:) .* ref_optical_image(:));
den3 = sqrt(sum(baseline_img(:).^2) * sum(ref_optical_image(:).^2)) + 1e-12;
ncc_base_vs_optical = num3 / den3;

dynamic_range_optical  = (max(ref_optical_image(:)) - min(ref_optical_image(:))) + 1e-12;
dynamic_range_baseline = (max(baseline_img(:))      - min(baseline_img(:)))      + 1e-12;

if exist('ssim','file') == 2
    ssim_opt_vs_optical  = ssim(optimized_img, ref_optical_image, 'DynamicRange', dynamic_range_optical);
    ssim_opt_vs_baseline = ssim(optimized_img, baseline_img,      'DynamicRange', dynamic_range_baseline);
    ssim_base_vs_optical = ssim(baseline_img,  ref_optical_image, 'DynamicRange', dynamic_range_optical);
else
    ssim_opt_vs_optical  = NaN;
    ssim_opt_vs_baseline = NaN;
    ssim_base_vs_optical = NaN;
end

if exist('psnr','file') == 2
    psnr_opt_vs_optical  = psnr(optimized_img, ref_optical_image, dynamic_range_optical);
    psnr_opt_vs_baseline = psnr(optimized_img, baseline_img,      dynamic_range_baseline);
    psnr_base_vs_optical = psnr(baseline_img,  ref_optical_image, dynamic_range_optical);
else
    psnr_opt_vs_optical  = 10 * log10((dynamic_range_optical^2)  / (mse_opt_vs_optical  + 1e-12));
    psnr_opt_vs_baseline = 10 * log10((dynamic_range_baseline^2) / (mse_opt_vs_baseline + 1e-12));
    psnr_base_vs_optical = 10 * log10((dynamic_range_optical^2)  / (mse_base_vs_optical + 1e-12));
end

fprintf('\n================ Image quality summary (REFERENCE) ================\n');
fprintf('\nDI-Optimized SAR vs Optical Reference:\n');
fprintf('MSE(DI-opt,Optical)   : %.4e\n', mse_opt_vs_optical);
fprintf('NCC(DI-opt,Optical)   : %.4f\n', ncc_opt_vs_optical);
fprintf('PSNR(DI-opt,Optical)  : %.2f dB\n', psnr_opt_vs_optical);
fprintf('SSIM(DI-opt,Optical)  : %.4f\n', ssim_opt_vs_optical);

fprintf('\nOptimized SAR vs Baseline SAR:\n');
fprintf('MSE(DI-opt,Baseline)  : %.4e\n', mse_opt_vs_baseline);
fprintf('NCC(DI-opt,Baseline)  : %.4f\n', ncc_opt_vs_baseline);
fprintf('PSNR(DI-opt,Baseline) : %.2f dB\n', psnr_opt_vs_baseline);
fprintf('SSIM(DI-opt,Baseline) : %.4f\n', ssim_opt_vs_baseline);

fprintf('\nBaseline SAR vs Optical Reference:\n');
fprintf('MSE(Baseline,Optical)  : %.4e\n', mse_base_vs_optical);
fprintf('NCC(Baseline,Optical)  : %.4f\n', ncc_base_vs_optical);
fprintf('PSNR(Baseline,Optical) : %.2f dB\n', psnr_base_vs_optical);
fprintf('SSIM(Baseline,Optical) : %.4f\n', ssim_base_vs_optical);

%% ============================================================
% VISUALIZATION
% ============================================================

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.75 0.75]);

subplot(2,2,1);
imagesc(xAxis, yAxis, baseline_img);
axis image;
set(gca, 'YDir', 'normal', 'Color', 'w', 'FontSize', 11, 'LineWidth', 1);
colormap(gca,'jet');
colorbar;
title(sprintf('Baseline %s', upper(sar_algo)));
xlabel('X (mm)');
ylabel('Y (mm)');

subplot(2,2,2);
imagesc(xAxis, yAxis, ref_optical_image);
axis image;
set(gca, 'YDir', 'normal', 'Color', 'w', 'FontSize', 11, 'LineWidth', 1);
colormap(gca,'gray');
colorbar;
title('Optical reference');
xlabel('X (mm)');
ylabel('Y (mm)');

subplot(2,2,3);
imagesc(xAxis, yAxis, optimized_img);
axis image;
set(gca, 'YDir', 'normal', 'Color', 'w', 'FontSize', 11, 'LineWidth', 1);
colormap(gca,'jet');
colorbar;
title(sprintf('Optimized %s', upper(sar_algo)));
xlabel('X (mm)');
ylabel('Y (mm)');

subplot(2,2,4);
imagesc(xAxis, yAxis, optimized_img - baseline_img);
axis image;
set(gca, 'YDir', 'normal', 'Color', 'w', 'FontSize', 11, 'LineWidth', 1);
colormap(gca,'jet');
colorbar;
title('Difference: Optimized - Baseline');
xlabel('X (mm)');
ylabel('Y (mm)');

sgtitle(sprintf('%s Optical-Guided Phase Calibration', upper(sar_algo)), ...
    'FontSize', 13, 'FontWeight', 'bold');

%% ============================================================
% FUNCTIONS
% ============================================================

function [loss, grad, current_img] = loss_and_grad_phaseOnly_MFA_RMA(sarData, PhaseVar, params)

    current_img = reconstruct_phaseOnly_MFA_RMA(sarData, PhaseVar, params);
    current_img = current_img / params.global_scale;

    L_img = mean((current_img - params.ref_optical_image).^2, 'all');

    lambda_phi = params.lambda_phi;
    L_reg_phi = lambda_phi * mean(PhaseVar.^2, 'all');

    lambda_smooth = params.lambda_smooth;
    dP_x = PhaseVar(:, 2:end) - PhaseVar(:, 1:end-1);
    dP_y = PhaseVar(2:end, :) - PhaseVar(1:end-1, :);
    L_reg_smooth = lambda_smooth * (mean(dP_x.^2, 'all') + mean(dP_y.^2, 'all'));

    lambda_dev = params.lambda_dev;
    operator_ref = params.operator_ref;
    operator_cur = params.operator_abs .* exp(1i * (params.operator_phase + PhaseVar));
    L_reg_dev = lambda_dev * mean(abs(operator_cur - operator_ref).^2, 'all');

    loss = real(L_img + L_reg_phi + L_reg_smooth + L_reg_dev);
    grad = dlgradient(loss, PhaseVar);
end

function img_dl = reconstruct_phaseOnly_MFA_RMA(sarData, PhaseVar, params)
    operator_cur = params.operator_abs .* exp(1i * (params.operator_phase + PhaseVar));

    switch upper(params.sar_algo)
        case 'MFA'
            [~, ~, img_dl, ~] = dlMFA_custom(sarData, operator_cur, params);
        case 'RMA'
            [~, ~, img_dl, ~] = dlRMA_custom(sarData, operator_cur, params);
        otherwise
            error('Unsupported sar_algo in reconstruct_phaseOnly_MFA_RMA.');
    end
end

function matchedFilter = refMF(params)
    c = physconst('lightspeed');
    x = params.dx * (-(params.nFFTspace-1)/2 : (params.nFFTspace-1)/2) * 1e-3;
    y = (params.dy * (-(params.nFFTspace-1)/2 : (params.nFFTspace-1)/2) * 1e-3).';
    z0_m = params.z0 * 1e-3;
    k = 2 * pi * params.F0 / c;

    matchedFilter = exp(-1i * 2 * k * sqrt(bsxfun(@plus, x.^2, y.^2) + z0_m^2));
end

function [xRangeT, yRangeT, trueImage_abs, trueImage_complx] = dlMFA_custom(sarData, matchedFilter, params)
    if isa(sarData,'dlarray') && ~isa(matchedFilter,'dlarray')
        matchedFilter = dlarray(matchedFilter);
    end

    [yPointM, xPointM] = size(sarData);
    [yPointF, xPointF] = size(matchedFilter);

    if (xPointF > xPointM)
        pad_x_pre  = floor((xPointF - xPointM) / 2);
        pad_x_post = ceil((xPointF - xPointM) / 2);
        sarData = cat(2, zeros(yPointM, pad_x_pre, 'like', sarData), ...
                         sarData, ...
                      zeros(yPointM, pad_x_post, 'like', sarData));
    end

    if (yPointF > yPointM)
        pad_y_pre  = floor((yPointF - yPointM) / 2);
        pad_y_post = ceil((yPointF - yPointM) / 2);
        sarData = cat(1, zeros(pad_y_pre, size(sarData,2), 'like', sarData), ...
                         sarData, ...
                      zeros(pad_y_post, size(sarData,2), 'like', sarData));
    end

    sarDataFFT        = fft(fft(sarData, [], 2), [], 1);
    matchedFilterFFT  = fft(fft(matchedFilter, [], 2), [], 1);
    trueImage_shifted = ifft(ifft(sarDataFFT .* matchedFilterFFT, [], 2), [], 1);
    trueImage         = fftshift(trueImage_shifted);

    [J, I] = size(trueImage);
    xij = round(params.bbox(1:2) / params.dx - 0.5 + I/2);
    ykl = round(params.bbox(3:4) / params.dy - 0.5 + J/2);

    trueImage_cropped = trueImage(ykl(1):ykl(2), xij(1):xij(2));
    trueImage_complx  = trueImage_cropped;
    trueImage_abs     = dlarray(abs(trueImage_cropped), 'SS');

    xRangeT = params.bbox(1) + (0:size(trueImage_abs,2)-1) * params.dx;
    yRangeT = params.bbox(3) + (0:size(trueImage_abs,1)-1) * params.dy;
end

function phaseFactor = buildRMAKernel(params)
    nFFTspace = params.nFFTspace;
    z0_use = params.z0;   % meters
    dx = params.dx;
    dy = params.dy;
    F0 = params.F0;

    c = physconst('lightspeed');
    k = 2*pi*F0/c;

    wSx = 2*pi / (dx * 1e-3);
    wSy = 2*pi / (dy * 1e-3);

    kX = linspace(-(wSx/2), (wSx/2), nFFTspace);
    kY = (linspace(-(wSy/2), (wSy/2), nFFTspace)).';

    arg = (2*k).^2 - bsxfun(@plus, kX.^2, kY.^2);
    evanescent_mask = arg < 0;

    K = sqrt(complex(arg, 0));
    phaseFactor = exp(-1i * z0_use * K);
    phaseFactor(evanescent_mask) = 0;
    phaseFactor = fftshift(fftshift(phaseFactor, 1), 2);
end

function [xRangeT, yRangeT, trueImage_abs, trueImage_complx] = dlRMA_custom(sarData, phaseFactor, params)
    if isa(sarData,'dlarray') && ~isa(phaseFactor,'dlarray')
        phaseFactor = dlarray(phaseFactor);
    end

    [yPointM, xPointM] = size(sarData);
    [yPointF, xPointF] = size(phaseFactor);

    if (xPointF > xPointM)
        pad_x_pre  = floor((xPointF - xPointM) / 2);
        pad_x_post = ceil((xPointF - xPointM) / 2);
        sarData = cat(2, zeros(yPointM, pad_x_pre, 'like', sarData), ...
                         sarData, ...
                      zeros(yPointM, pad_x_post, 'like', sarData));
    elseif (xPointM > xPointF)
        pad_x_pre  = floor((xPointM - xPointF) / 2);
        pad_x_post = ceil((xPointM - xPointF) / 2);
        phaseFactor = cat(2, zeros(yPointF, pad_x_pre, 'like', phaseFactor), ...
                             phaseFactor, ...
                          zeros(yPointF, pad_x_post, 'like', phaseFactor));
    end

    if (yPointF > yPointM)
        pad_y_pre  = floor((yPointF - yPointM) / 2);
        pad_y_post = ceil((yPointF - yPointM) / 2);
        sarData = cat(1, zeros(pad_y_pre, size(sarData,2), 'like', sarData), ...
                         sarData, ...
                      zeros(pad_y_post, size(sarData,2), 'like', sarData));
    elseif (yPointM > yPointF)
        pad_y_pre  = floor((yPointM - yPointF) / 2);
        pad_y_post = ceil((yPointM - yPointF) / 2);
        phaseFactor = cat(1, zeros(pad_y_pre, size(phaseFactor,2), 'like', phaseFactor), ...
                             phaseFactor, ...
                          zeros(pad_y_post, size(phaseFactor,2), 'like', phaseFactor));
    end

    sarDataFFT = fft(fft(sarData, [], 2), [], 1);
    trueImage  = ifft(ifft(sarDataFFT .* phaseFactor, [], 2), [], 1);

    [J, I] = size(trueImage);
    xij = round(params.bbox(1:2) / params.dx - 0.5 + I/2);
    ykl = round(params.bbox(3:4) / params.dy - 0.5 + J/2);

    trueImage_cropped = trueImage(ykl(1):ykl(2), xij(1):xij(2));
    trueImage_complx  = trueImage_cropped;
    trueImage_abs     = dlarray(abs(trueImage_cropped), 'SS');

    xRangeT = params.bbox(1) + (0:size(trueImage_abs,2)-1) * params.dx;
    yRangeT = params.bbox(3) + (0:size(trueImage_abs,1)-1) * params.dy;
end

function M = image_quality_metrics(img)
    A = abs(img);
    A2 = A.^2;

    M.peak = max(A(:));

    denom2 = sum(A2(:))^2 + eps;
    M.contrast = sqrt(numel(A) * sum(A(:).^4) / denom2);

    p = A2 / (sum(A2(:)) + eps);
    M.entropy = -sum(p(:) .* log(p(:) + eps));
end