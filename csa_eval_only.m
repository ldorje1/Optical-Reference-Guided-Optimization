clc; clear; close all;

%% ========================================================================
% Conventional CSA for Near-Field SAR/mmWave Imaging
% Uses SAME dataset convention/properties as main optical-reference code
% No optical optimization. Evaluation only.
%% ========================================================================

%% -------------------- USER SELECTION --------------------
sar_algo = 'CSA';
dataset  = 3;   % 1..6
% Dataset convention:
%   1 -> YANIK
%   2 -> 3DRIED case 1
%   3 -> 3DRIED case 2
%   4 -> OURS case 1
%   5 -> OURS case 2
%   6 -> OURS screw_driver

raw_dataDir    = fullfile(pwd, 'sar_raw_data');
optical_ref_dir = fullfile(pwd, 'optical_refs');

% Image grid for CSA reconstruction
A = 80;   % image width
B = 80;   % image height

% FFT size
nFFTtime = 1024;

% Constants
c0 = physconst('lightspeed');

%% ========================================================================
% DATASET-SPECIFIC PREPROCESSING
%% ========================================================================

if dataset == 1
    % ==============================================================
    % DATASET 1: YANIK
    % ==============================================================
    obj_name = "yanik";
    dataName = "rawData3D_simple2D.mat";

    bbox = [-70 80 -80 70];
    z0   = 280;              % mm

    nFFTtime  = 1024;
    nFFTspace = 1024;

    F0 = 77e9;
    FS = 9121e3;
    K0 = 63.343e12;
    tI = 4.5225e-10;

    optical_ref_file    = 'ref_optical_image_yanik.mat';
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
    % DATASET 2 or 3: 3DRIED
    % ==============================================================
    obj_name = "3dried";

    if dataset == 2
        dataName = "rawdata_gun_3dried.mat";
        bbox = [-150 250 -320 80];
        optical_ref_file = 'ref_optical_image_gun_3dried.mat';
    else
        dataName = "rawdata_knife_3dried.mat";
        bbox = [-160 250 -210 260];
        optical_ref_file = 'ref_optical_image_knife_3dried.mat';
    end

    optical_ref_varname = 'ref_optical_image';
    invert_optical_reference = true;

    F0 = (77 + 1.8) * 1e9;
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
        Sr0 = Raw_echo(1:Num:end, :);
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

    sarRawData = Echo;     % Nx x Nz x Nsamp
    Nsamp      = num_sample;
    M          = Nz;
    N          = Nx;

    ID_select = 17;
    z0 = (c0/2) * ((ID_select/(K0*(1/FS)*Nsamp)) - tI) * 1000;   % mm

elseif dataset == 4 || dataset == 5 || dataset == 6
    % ==============================================================
    % DATASET 4, 5, or 6: OURS
    % ==============================================================
    obj_name = "ours";

    if dataset == 4
        dataName = "rawdata_knife_ours.mat";
        dx = 1; dy = 1;
        z0 = 185;
        bbox = [-200 200 -200 200];
        optical_ref_file = 'ref_optical_image_knife_ours.mat';
    elseif dataset == 5
        dataName = "rawdata_plier_ours.mat";
        dx = 1; dy = 2;
        z0 = 210;
        bbox = [-200 200 -200 200];
        optical_ref_file = 'ref_optical_image_plier_ours.mat';
    else
        dataName = "rawdata_screw_driver_ours.mat";
        dx = 1; dy = 2;
        z0 = 230;
        bbox = [-200 200 -200 200];
        optical_ref_file = 'ref_optical_image_screw_driver_ours.mat';
    end

    optical_ref_varname = 'ref_optical_image';
    invert_optical_reference = true;

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

    [Nsamp, M, N] = size(sarRawData);

else
    error('Unknown dataset selection. Choose dataset = 1, 2, 3, 4, 5, or 6.');
end

fprintf('Selected dataset   : %d\n', dataset);
fprintf('Selected family    : %s\n', upper(obj_name));
fprintf('Selected algorithm : %s\n', upper(sar_algo));
fprintf('dx = %.3f mm, dy = %.3f mm, z0 = %.3f mm, FS = %.3f MHz\n', ...
    dx, dy, z0, FS/1e6);

%% ========================================================================
% LOAD OPTICAL REFERENCE
%% ========================================================================

opticalStruct = load(fullfile(optical_ref_dir, optical_ref_file));
if ~isfield(opticalStruct, optical_ref_varname)
    error('Variable "%s" not found in file "%s".', ...
        optical_ref_varname, optical_ref_file);
end

ref_optical_image = double(opticalStruct.(optical_ref_varname));

if invert_optical_reference
    ref_optical_image = 1 - ref_optical_image;
end

%% ========================================================================
% RANGE GATING
%% ========================================================================


if dataset == 2 || dataset == 3
    % 3DRIED path follows same style as main code
    rawDataFFT = fft(sarRawData, Nsamp, 3);         % Nx x Nz x Nsamp
    k0_range_bin = 17;
    sarData = squeeze(rawDataFFT(:, :, k0_range_bin)).';   % Nz x Nx

    for ii = 2:2:size(sarData,1)
        sarData(ii,:) = fliplr(sarData(ii,:));
    end

else
    % YANIK and OURS path
    Echo = permute(sarRawData, [3, 2, 1]);          % Nx x Ny x Nsamp
    num_sample = size(Echo, 3);

    rawDataFFT = fft(Echo, num_sample, 3);

    k0_range_bin = round(K0 * (1 / FS) * (2 * z0 * 1e-3 / c0 + tI) * num_sample);
    
    %E = squeeze(sum(sum(abs(rawDataFFT).^2, 1), 2));
    %[~, k0_range_bin] = max(E);

    sarData = squeeze(rawDataFFT(:, :, k0_range_bin + 1)).';

    if dataset == 4 || dataset == 5 || dataset == 6 
        for ii = 2:2:size(sarData,1)
            sarData(ii,:) = fliplr(sarData(ii,:));
        end
    end
end
    


%% ========================================================================
% PARAMS
%% ========================================================================

params = struct;
params.z0           = z0;
params.dx           = dx;
params.dy           = dy;
params.bbox         = bbox;
params.Nsamp        = Nsamp;
params.nFFTtime     = nFFTtime;
params.N            = N;
params.M            = M;
params.A            = A;
params.B            = B;
params.F0           = F0;
params.K0           = K0;
params.FS           = FS;
params.tI           = tI;
params.k0_range_bin = k0_range_bin;
params.sar_algo     = upper(sar_algo);

% CSA / SBRIM parameters
params.lambda0_csa  = 1e-4;
params.p_csa        = 1.0;
params.eta_csa      = 1e-5;
params.maxIter_csa  = 50;
params.epsilon0_csa = 1e-4;

%% ========================================================================
% BUILD CSA MODEL
%% ========================================================================

H_csa        = dlCSA_H_matrix(params);
params.H_csa = H_csa;

%% ========================================================================
% RECONSTRUCT
%% ========================================================================

[xAxis, yAxis, img_dl, img_cmplx_dl, alpha_hat_dl] = dlCSA_fast(sarData, params);

img_eval  = abs(gather(extractdata(img_dl)));
img_eval = fliplr(img_eval);
img_cmplx = gather(extractdata(img_cmplx_dl));
img_plot  = img_eval / (max(img_eval(:)) + eps);

%% ========================================================================
% MATCH OPTICAL SIZE + NORMALIZATION
%% ========================================================================

if ~isequal(size(ref_optical_image), size(img_eval))
    fprintf('Resizing optical reference %s -> %s\n', ...
        mat2str(size(ref_optical_image)), mat2str(size(img_eval)));

    if exist('imresize','file') == 2
        ref_optical_image = imresize(ref_optical_image, size(img_eval), 'bilinear');
    else
        [X, Y]  = meshgrid(linspace(0,1,size(ref_optical_image,2)), ...
                           linspace(0,1,size(ref_optical_image,1)));
        [Xq,Yq] = meshgrid(linspace(0,1,size(img_eval,2)), ...
                           linspace(0,1,size(img_eval,1)));
        ref_optical_image = interp2(X, Y, ref_optical_image, Xq, Yq, 'linear', 0);
    end
end

global_scale = max(abs(img_eval(:))) + 1e-12;
img_eval_norm = img_eval / global_scale;
ref_optical_image = ref_optical_image / (max(abs(ref_optical_image(:))) + 1e-12);

%% ========================================================================
% SHOW CSA + OPTICAL
%% ========================================================================

figure;
subplot(1,2,1);
imagesc(xAxis, yAxis, img_eval_norm);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca,'jet');
colorbar;
title('CSA Reconstruction');

subplot(1,2,2);
imagesc(xAxis, yAxis, ref_optical_image);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca,'gray');
colorbar;
title('Optical Reference');

%% ========================================================================
% NO-REFERENCE IMAGE QUALITY METRICS
%% ========================================================================

Q = image_quality_metrics(img_eval_norm);

fprintf('\n================ Image quality summary (NO-REFERENCE) ================\n');
fprintf('CSA:\n');
fprintf('  Contrast : %.6f\n', Q.contrast);
fprintf('  Entropy  : %.6f\n', Q.entropy);
fprintf('  Peak     : %.6e\n', Q.peak);
fprintf('========================================================================\n');

%% ========================================================================
% REFERENCE-BASED METRICS: CSA vs OPTICAL REFERENCE
% Same metric-only normalization as DI-MFA/DI-RMA code
%% ========================================================================

csa_m = img_eval_norm;
ref_m = ref_optical_image;

% Metric-only min-max normalization
% This does not change displayed images.
csa_m = csa_m - min(csa_m(:));
csa_m = csa_m / (max(csa_m(:)) + 1e-12);

ref_m = ref_m - min(ref_m(:));
ref_m = ref_m / (max(ref_m(:)) + 1e-12);

% ============================================================
% TBR METRIC USING OPTICAL-REFERENCE TARGET MASK
% ============================================================

target_mask = create_target_mask_from_optical(ref_m);

tbr_csa = compute_tbr_from_mask(csa_m, target_mask);

mse_csa_vs_optical = mean((csa_m(:) - ref_m(:)).^2);

num = sum(csa_m(:) .* ref_m(:));
den = sqrt(sum(csa_m(:).^2) * sum(ref_m(:).^2)) + 1e-12;
cos_csa_vs_optical = num / den;

if exist('ssim','file') == 2
    ssim_csa_vs_optical = ssim(csa_m, ref_m, 'DynamicRange', 1);
else
    ssim_csa_vs_optical = NaN;
end

if exist('psnr','file') == 2
    psnr_csa_vs_optical = psnr(csa_m, ref_m, 1);
else
    psnr_csa_vs_optical = 10 * log10(1 / (mse_csa_vs_optical + 1e-12));
end

fprintf('\n================ Image quality summary (REFERENCE) ================\n');
fprintf('\nCSA vs Optical Reference:\n');
fprintf('MSE(CSA,Optical)      : %.4e\n', mse_csa_vs_optical);
fprintf('CosSim(CSA,Optical)   : %.4f\n', cos_csa_vs_optical);
fprintf('PSNR(CSA,Optical)     : %.2f dB\n', psnr_csa_vs_optical);
fprintf('SSIM(CSA,Optical)     : %.4f\n', ssim_csa_vs_optical);
fprintf('TBR(CSA)              : %.2f dB\n', tbr_csa);

%% ========================================================================
% VISUALIZATION
%% ========================================================================

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.75 0.75]);
xlims = [xAxis(1) xAxis(end)];
ylims = [yAxis(1) yAxis(end)];

subplot(1,2,1);
imagesc(xAxis, yAxis, img_eval_norm);
set(gca,'YDir','normal');
axis image;
xlim(xlims); ylim(ylims);
colormap(gca,'jet'); colorbar;
title('CSA Reconstruction');

subplot(1,2,2);
imagesc(xAxis, yAxis, ref_optical_image);
set(gca,'YDir','normal');
axis image;
xlim(xlims); ylim(ylims);
colormap(gca,'gray'); colorbar;
title('Optical Reference');

sgtitle(sprintf('CSA Evaluation | Dataset %d', dataset), ...
    'FontSize', 13, 'FontWeight', 'bold');

%% ========================================================================
% FUNCTIONS
%% ========================================================================

function [xRangeT, yRangeT, trueImage_abs, trueImage_complx, alpha_hat_dl] = dlCSA_fast(sarData, params)

    if isa(sarData, 'dlarray')
        sarData_num = double(extractdata(sarData));
    else
        sarData_num = double(sarData);
    end

    if isa(params.H_csa, 'dlarray')
        H_num = double(extractdata(params.H_csa));
    else
        H_num = double(params.H_csa);
    end

    ys_num = sarData_num(:);

    alpha_hat_num = CSA_SBRIM_numeric(ys_num, H_num, ...
                                      params.lambda0_csa, ...
                                      params.p_csa, ...
                                      params.eta_csa, ...
                                      params.maxIter_csa, ...
                                      params.epsilon0_csa);

    B = params.B;
    A = params.A;
    alpha_img_num = reshape(alpha_hat_num, B, A);

    trueImage_complx = dlarray(alpha_img_num);
    trueImage_abs    = dlarray(abs(alpha_img_num), 'SS');
    alpha_hat_dl     = dlarray(alpha_hat_num);

    xRangeT = linspace(params.bbox(1), params.bbox(2), A);
    yRangeT = linspace(params.bbox(3), params.bbox(4), B);
end

function alpha_hat = CSA_SBRIM_numeric(ys, H, lambda0, p, eta, maxIter, epsilon0)

    ys = double(ys);
    H  = double(H);

    [M_meas, ~] = size(H);

    HtH  = H' * H;
    Htys = H' * ys;

    alpha_hat = Htys;
    r = Inf;
    n = 0;
    beta_n = 1;

    fprintf('Starting SBRIM (numeric), p=%.2f...\n', p);

    while (r >= epsilon0) && (n < maxIter)
        n = n + 1;
        alpha_prev = alpha_hat;

        alpha_sq_plus_eta = abs(alpha_prev).^2 + eta;
        lambda_diag       = (p / 2) * (alpha_sq_plus_eta).^(p/2 - 1);

        A_mat = HtH + lambda0 * beta_n * diag(lambda_diag);
        alpha_hat = A_mat \ Htys;

        residual = ys - H * alpha_hat;
        beta_n   = sum(abs(residual).^2) / M_meas;

        denom = norm(alpha_hat);
        if denom < eps
            r = 0;
        else
            r = norm(alpha_hat - alpha_prev) / denom;
        end

        if mod(n, 10) == 0 || n == 1
            fprintf('Iter %d: r=%.4e, beta=%.4e\n', n, r, beta_n);
        end
    end

    if n == maxIter
        fprintf('Warning: SBRIM reached maxIter=%d (r=%.4e)\n', maxIter, r);
    else
        fprintf('SBRIM converged in %d iters (r=%.4e)\n', n, r);
    end
end

function H = dlCSA_H_matrix(params)

    Ny = params.M;
    Nx = params.N;
    A  = params.A;
    B  = params.B;

    c0    = physconst('lightspeed');
    F0    = params.F0;
    z0_mm = params.z0;
    dx    = params.dx;
    dy    = params.dy;
    bbox  = params.bbox;

    z0_m   = z0_mm * 1e-3;
    dxm    = dx   * 1e-3;
    dym    = dy   * 1e-3;
    bbox_m = bbox * 1e-3;

    k   = 2 * pi * F0 / c0;
    cst = 1i * 2 * k;
    z2  = z0_m^2;

    wh1 = linspace(bbox_m(1), bbox_m(2), A);
    wh2 = linspace(bbox_m(3), bbox_m(4), B);

    NM = Ny * Nx;
    BA = A * B;

    H_val = complex(zeros(NM, BA));

    fprintf('Building CSA H matrix (%d x %d)... ', NM, BA);
    tic;
    for i = 1:NM
        iy = mod(i-1, Ny);
        ix = (i-1-iy) / Ny;

        sx_i = (ix + 0.5 - Nx/2) * dxm;
        sy_i = (iy + 0.5 - Ny/2) * dym;

        for j = 1:BA
            jy = mod(j-1, B);
            jx = (j-1-jy) / B;

            px = wh1(jx+1);
            py = wh2(jy+1);

            dist2      = (sx_i - px)^2 + (sy_i - py)^2 + z2;
            H_val(i,j) = exp(cst * sqrt(dist2));
        end
    end
    fprintf('%.3f sec\n', toc);

    H = dlarray(H_val);
end

function M = image_quality_metrics(img)
    A = abs(double(img));
    A2 = A.^2;
    M.peak = max(A(:));
    denom2 = sum(A2(:))^2 + eps;
    M.contrast = sqrt(max(numel(A) * sum(A(:).^4) / denom2 - 1, 0));
    p = A2 / (sum(A2(:)) + eps);
    M.entropy = -sum(p(:) .* log(p(:) + eps));
end


function target_mask = create_target_mask_from_optical(ref_m)

    ref_m = double(ref_m);

    ref_m = ref_m - min(ref_m(:));
    ref_m = ref_m / (max(ref_m(:)) + eps);

    threshold = 0.30;
    target_mask = ref_m > threshold;

    if exist('imfill','file') == 2
        target_mask = imfill(target_mask, 'holes');
    end

    if exist('bwareaopen','file') == 2
        target_mask = bwareaopen(target_mask, 20);
    end

    if nnz(target_mask) < 10
        if exist('graythresh','file') == 2
            threshold = graythresh(ref_m);
            target_mask = ref_m > threshold;

            if exist('imfill','file') == 2
                target_mask = imfill(target_mask, 'holes');
            end

            if exist('bwareaopen','file') == 2
                target_mask = bwareaopen(target_mask, 20);
            end
        end
    end

    if nnz(target_mask) < 10
        warning('Target mask is very small. TBR may be unreliable.');
    end
end

function TBR_dB = compute_tbr_from_mask(img, target_mask)

    I = abs(double(img));

    I = I - min(I(:));
    I = I / (max(I(:)) + eps);

    target_mask = logical(target_mask);
    background_mask = ~target_mask;

    if nnz(target_mask) == 0 || nnz(background_mask) == 0
        TBR_dB = NaN;
        warning('Invalid target/background mask. TBR set to NaN.');
        return;
    end

    Pt = mean(I(target_mask).^2);
    Pb = mean(I(background_mask).^2);

    TBR_dB = 10 * log10((Pt + eps) / (Pb + eps));
end