clc; clear; close all;

%% ========================================================================
% Conventional BPA / Enhanced BPA / LIA for Near-Field SAR/mmWave Imaging
% Uses SAME dataset convention/properties as main optical-reference code
% Includes baseline-vs-optical evaluation only
% No optimization.
%% ========================================================================

%% -------------------- USER SELECTION --------------------
sar_algo = 'ENHANCEDBPA';   % 'BPA', 'ENHANCEDBPA', or 'LIA'
dataset  = 6;       % 1..6

% Dataset convention:
%   1 -> YANIK
%   2 -> 3DRIED case 1
%   3 -> 3DRIED case 2
%   4 -> OURS case 1
%   5 -> OURS case 2
%   6 -> OURS screw_driver

raw_dataDir    = fullfile(pwd, 'sar_raw_data');
optical_ref_dir = fullfile(pwd, 'optical_refs');

% Imaging box / grid
switch upper(sar_algo)
    case 'LIA'
        A = 60;
        B = 60;

    case {'BPA','ENHANCEDBPA'}
        A = 80;
        B = 80;

    otherwise
        error('sar_algo must be BPA, LIA, or EnhancedBPA.');
end

% FFT sizes
nFFTtime = 1024;

% Radar constants
c0 = physconst('lightspeed');

% Enhanced BPA parameters
params_enh.beam_Mx     = 7;
params_enh.beam_My     = 7;
params_enh.beam_dx     = [];   % leave empty -> use dx
params_enh.beam_dy     = [];   % leave empty -> use dy
params_enh.use_window  = true;
params_enh.window_type = 'none';   % 'hann' | 'hamming' | 'none'

%% ========================================================================
% DATASET-SPECIFIC PREPROCESSING (matched to main code)
%% ========================================================================

if dataset == 1
    % ==============================================================
    % DATASET 1 YANIK
    % ==============================================================
    obj_name = "yanik";
    dataName = "rawData3D_simple2D.mat";

    bbox = [-70 80 -80 70];
    z0   = 280;                        % mm

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
    % DATASET 2 or 3 (3DRIED)
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
    M = Nz;
    N = Nx;

    ID_select = 17;
    z0 = (c0/2) * ((ID_select/(K0*(1/FS)*Nsamp)) - tI) * 1000;   % mm

elseif dataset == 4 || dataset == 5 || dataset == 6
    % ==============================================================
    % DATASET 4, 5, or 6 OURS
    % ==============================================================
    obj_name = "ours";

    if dataset == 4
        dataName = "rawdata_knife_ours.mat";
        optical_ref_file = 'ref_optical_image_knife_ours.mat';
        dx = 1; dy = 1;
        z0   = 185;
        bbox = [-200 200 -200 200];
    elseif dataset == 5
        dataName = "rawdata_plier_ours.mat";
        optical_ref_file = 'ref_optical_image_plier_ours.mat';
        dx = 1; dy = 2;
        z0   = 210;
        bbox = [-200 200 -200 200];
    else
        dataName = "rawdata_screw_driver_ours.mat";
        optical_ref_file = 'ref_optical_image_screw_driver_ours.mat';
        dx = 1; dy = 2;
        z0   = 230;
        bbox = [-200 200 -200 200];
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

if isempty(params_enh.beam_dx), params_enh.beam_dx = dx; end
if isempty(params_enh.beam_dy), params_enh.beam_dy = dy; end

fprintf('Selected dataset   : %d\n', dataset);
fprintf('Selected family    : %s\n', upper(obj_name));
fprintf('Selected algorithm : %s\n', upper(sar_algo));
fprintf('dx = %.4f mm, dy = %.4f mm, z0 = %.4f mm, FS = %.4f MHz\n', ...
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
    rawDataFFT = fft(sarRawData, Nsamp, 3);         % Nx x Nz x Nsamp
    k0_range_bin = 17;
    sarData = squeeze(rawDataFFT(:, :, k0_range_bin)).';   % Nz x Nx

    for ii = 2:2:size(sarData,1)
        sarData(ii,:) = fliplr(sarData(ii,:));
    end

else
    rawDataFFT   = fft(sarRawData, nFFTtime, 1);
    k0_range_bin = round(K0 * (1 / FS) * (2 * z0 * 1e-3 / c0 + tI) * nFFTtime);
    sarData      = squeeze(rawDataFFT(k0_range_bin + 1, :, :));

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
params.A_bpa        = A;
params.B_bpa        = B;
params.A_img        = A;
params.B_img        = B;
params.F0           = F0;
params.K0           = K0;
params.FS           = FS;
params.tI           = tI;
params.k0_range_bin = k0_range_bin;
params.sar_algo     = upper(sar_algo);
params.beam_Mx      = params_enh.beam_Mx;
params.beam_My      = params_enh.beam_My;
params.beam_dx      = params_enh.beam_dx;
params.beam_dy      = params_enh.beam_dy;
params.use_window   = params_enh.use_window;
params.window_type  = params_enh.window_type;

%% ========================================================================
% RECONSTRUCT
%% ========================================================================

switch upper(sar_algo)
    case 'BPA'
        H_bpa = dlBPA_H_matrix(params);
        params.H_bpa = H_bpa;
        [xAxis, yAxis, img_dl, img_cmplx_dl] = dlBPA(sarData, params, H_bpa);
        method_name = 'Conventional BPA';

    case 'LIA'
        H_bpa = dlBPA_H_matrix(params);
        params.H_bpa = H_bpa;
        NM = M * N;
        kk = min(40000, NM);
        rng(1000);
        params.py = sort(randperm(NM, kk));
        [xAxis, yAxis, img_dl, img_cmplx_dl] = dlLIA(sarData, params, H_bpa);
        method_name = 'Conventional LIA';

    case 'ENHANCEDBPA'
        [xAxis, yAxis, img_abs, img_cmplx] = enhancedBPA2D(sarData, params);
        img_dl = dlarray(img_abs, 'SS');
        img_cmplx_dl = dlarray(img_cmplx);
        method_name = 'Enhanced BPA';

    otherwise
        error('sar_algo must be BPA, EnhancedBPA, or LIA.');
end

baseline_img = abs(gather(extractdata(img_dl)));
img_cmplx    = gather(extractdata(img_cmplx_dl));

%% ========================================================================
% MATCH OPTICAL SIZE + GLOBAL NORMALIZATION
%% ========================================================================

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
baseline_img = baseline_img / global_scale;
ref_optical_image = ref_optical_image / (max(abs(ref_optical_image(:))) + 1e-12);

%% ========================================================================
% SHOW BASELINE + OPTICAL
%% ========================================================================

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

%% ========================================================================
% NO-REFERENCE IMAGE QUALITY METRICS
%% ========================================================================

Q = image_quality_metrics(baseline_img);

fprintf('\n================ Image quality summary (NO-REFERENCE) ================\n');
fprintf('Baseline %s:\n', upper(sar_algo));
fprintf('  Contrast : %.6f\n', Q.contrast);
fprintf('  Entropy  : %.6f\n', Q.entropy);
fprintf('  Peak     : %.6e\n', Q.peak);
fprintf('========================================================================\n');

%% ========================================================================
% REFERENCE-BASED METRICS: BASELINE SAR vs OPTICAL REFERENCE
% Same metric-only normalization as DI-MFA/DI-RMA and CSA code
%% ========================================================================

base_m = baseline_img;
ref_m  = ref_optical_image;

% Metric-only min-max normalization
% This does not change displayed images.
base_m = base_m - min(base_m(:));
base_m = base_m / (max(base_m(:)) + 1e-12);

ref_m = ref_m - min(ref_m(:));
ref_m = ref_m / (max(ref_m(:)) + 1e-12);

% ============================================================
% TBR METRIC USING OPTICAL-REFERENCE TARGET MASK
% ============================================================

target_mask = create_target_mask_from_optical(ref_m);

tbr_base = compute_tbr_from_mask(base_m, target_mask);

mse_base_vs_optical = mean((base_m(:) - ref_m(:)).^2);

num = sum(base_m(:) .* ref_m(:));
den = sqrt(sum(base_m(:).^2) * sum(ref_m(:).^2)) + 1e-12;
cos_base_vs_optical = num / den;

if exist('ssim','file') == 2
    ssim_base_vs_optical = ssim(base_m, ref_m, 'DynamicRange', 1);
else
    ssim_base_vs_optical = NaN;
end

if exist('psnr','file') == 2
    psnr_base_vs_optical = psnr(base_m, ref_m, 1);
else
    psnr_base_vs_optical = 10 * log10(1 / (mse_base_vs_optical + 1e-12));
end

fprintf('\n================ Image quality summary (REFERENCE) ================\n');
fprintf('\nBaseline SAR vs Optical Reference:\n');
fprintf('MSE(Baseline,Optical)     : %.4e\n', mse_base_vs_optical);
fprintf('CosSim(Baseline,Optical)  : %.4f\n', cos_base_vs_optical);
fprintf('PSNR(Baseline,Optical)    : %.2f dB\n', psnr_base_vs_optical);
fprintf('SSIM(Baseline,Optical)    : %.4f\n', ssim_base_vs_optical);
fprintf('TBR(Baseline)             : %.2f dB\n', tbr_base);

%% ========================================================================
% VISUALIZATION
%% ========================================================================

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.75 0.75]);
xlims = [xAxis(1) xAxis(end)];
ylims = [yAxis(1) yAxis(end)];

subplot(1,2,1);
imagesc(xAxis, yAxis, baseline_img);
set(gca,'YDir','normal');
axis image;
xlim(xlims); ylim(ylims);
colormap(gca,'jet'); colorbar;
title(sprintf('Baseline %s', upper(sar_algo)));

subplot(1,2,2);
imagesc(xAxis, yAxis, ref_optical_image);
set(gca,'YDir','normal');
axis image;
xlim(xlims); ylim(ylims);
colormap(gca,'gray'); colorbar;
title('Optical reference');

sgtitle(sprintf('%s Evaluation | Dataset %d', method_name, dataset), ...
    'FontSize', 13, 'FontWeight', 'bold');

%% ========================================================================
% FUNCTIONS
%% ========================================================================

%%%%%%%%%%%%%%%%%%% BPA
function [xRangeT, yRangeT, trueImage_abs, trueImage_complx] = dlBPA(sarData, params, H)
    A = params.A_bpa;
    B = params.B_bpa;

    if ~isa(sarData, 'dlarray')
        sarData = dlarray(sarData);
    end

    y = reshape(sarData, [], 1);
    xd = H' * y;
    xdi = reshape(xd, B, A);
    trueImage_cropped = fliplr(xdi);
    trueImage_complx = trueImage_cropped;
    trueImage_abs = dlarray(abs(trueImage_complx), 'SS');

    xRangeT = linspace(params.bbox(1), params.bbox(2), A);
    yRangeT = linspace(params.bbox(3), params.bbox(4), B);
end

function H = dlBPA_H_matrix(params)
    A = params.A_bpa;
    B = params.B_bpa;

    c0 = physconst('lightspeed');
    F0 = params.F0;
    z0_mm = params.z0;
    dx = params.dx;
    dy = params.dy;
    bbox = params.bbox;

    z0_m = z0_mm * 1e-3;
    dxm = dx * 1e-3;
    dym = dy * 1e-3;
    bbox_m = bbox * 1e-3;

    k   = 2*pi*F0/c0;
    cst = 1i * 2 * k;
    z2  = z0_m^2;

    wh1 = linspace(bbox_m(1), bbox_m(2), A);
    wh2 = linspace(bbox_m(3), bbox_m(4), B);

    Ny = params.M;
    Nx = params.N;
    NM = Ny * Nx;
    BA = A * B;

    H_val = complex(zeros(NM, BA));
    fprintf('Building H matrix (%d x %d)... ', NM, BA);
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

            dist2 = (sx_i - px)^2 + (sy_i - py)^2 + z2;
            H_val(i,j) = exp(cst * sqrt(dist2));
        end
    end

    fprintf('%.3f sec\n', toc);
    H = dlarray(H_val);
end

%%%%%%%%%%%%%%%%%%% LIA
function [xRangeT, yRangeT, trueImage_abs, trueImage_complx] = dlLIA(sarData, params, H_bpa)
    if ~isa(sarData, 'dlarray')
        sarData = dlarray(sarData);
    end
    if ~isa(H_bpa, 'dlarray')
        H_bpa = dlarray(H_bpa);
    end

    M  = params.M;
    N  = params.N;
    A  = params.A_bpa;
    B  = params.B_bpa;
    py = params.py;

    rd_full = reshape(sarData, [], 1);
    rd      = rd_full(py);

    Hp = H_bpa(py, :);
    BA = A * B;

    di = 0.01;
    G  = di * (Hp' * Hp);
    xd = di * (Hp' * rd);

    for j = 1:BA
        Gj    = G(:, j);
        denom = 1 + G(j, j);
        temp  = Gj / denom;
        xd = xd - temp * xd(j);
        G  = G  - temp * G(j, :);
    end

    diagG = G(1:BA+1:BA*BA);
    diagG = reshape(diagG, [BA, 1]);
    xd = xd ./ diagG;

    xdi = fliplr(reshape(xd, B, A));
    trueImage_complx = xdi;
    trueImage_abs = dlarray(abs(trueImage_complx), 'SS');

    xRangeT = linspace(params.bbox(1), params.bbox(2), A);
    yRangeT = linspace(params.bbox(3), params.bbox(4), B);
end


%%%%%%%%%%%%%%%%%%% BPA
function [xRange, yRange, img_abs, img_cplx] = enhancedBPA2D(sarData, params)
    c0 = 299792458;
    lambda = c0 / params.F0;
    k = 2*pi / lambda;

    [X, Y, xRange, yRange] = make_image_grid(params);
    [sx, sy, winVec] = make_sensor_grid_and_window(params, size(sarData,1), size(sarData,2));

    z0_m = params.z0 * 1e-3;
    sarVec = reshape(sarData, [], 1);

    img = complex(zeros(size(X)));

    for n = 1:numel(sarVec)
        dx_rel = X - sx(n);
        dy_rel = Y - sy(n);
        R = sqrt(dx_rel.^2 + dy_rel.^2 + z0_m^2);

        phaseTerm = exp(-1i * 2 * k * R);
        G = beam_weight_map(dx_rel, dy_rel, R, params, lambda);

        img = img + winVec(n) * sarVec(n) .* phaseTerm .* G;
    end

    img_cplx = fliplr(img);
    img_abs  = abs(img_cplx);
end

function G = beam_weight_map(dx_rel, dy_rel, R, params, lambda)
    sin_theta_x = dx_rel ./ (R + eps);
    sin_theta_y = dy_rel ./ (R + eps);

    Gx = normalized_array_factor(sin_theta_x, params.beam_Mx, params.beam_dx * 1e-3, lambda);
    Gy = normalized_array_factor(sin_theta_y, params.beam_My, params.beam_dy * 1e-3, lambda);

    G = Gx .* Gy;
    G = G / (max(G(:)) + eps);
end

function G = normalized_array_factor(sinTheta, M_eff, d_elem, lambda)
    if M_eff <= 1
        G = ones(size(sinTheta));
        return;
    end

    beta = 2*pi * d_elem / lambda .* sinTheta;
    num = sin(M_eff * beta / 2);
    den = M_eff * sin(beta / 2);

    G = abs(num ./ (den + eps));
    smallMask = abs(beta) < 1e-10;
    G(smallMask) = 1;
    G = min(G, 1);
end

function [X, Y, xRange, yRange] = make_image_grid(params)
    xRange = linspace(params.bbox(1), params.bbox(2), params.A_img);
    yRange = linspace(params.bbox(3), params.bbox(4), params.B_img);
    [Xmm, Ymm] = meshgrid(xRange, yRange);
    X = Xmm * 1e-3;
    Y = Ymm * 1e-3;
end

function [sx, sy, winVec] = make_sensor_grid_and_window(params, Ny_ap, Nx_ap)
    dxm = params.dx * 1e-3;
    dym = params.dy * 1e-3;

    [colGrid, rowGrid] = meshgrid(0:Nx_ap-1, 0:Ny_ap-1);
    sx = ((colGrid(:) + 0.5) - Nx_ap/2) * dxm;
    sy = ((rowGrid(:) + 0.5) - Ny_ap/2) * dym;

    if params.use_window
        wx = local_window(Nx_ap, params.window_type);
        wy = local_window(Ny_ap, params.window_type);
        [WX, WY] = meshgrid(wx, wy);
        W = WX .* WY;
    else
        W = ones(Ny_ap, Nx_ap);
    end

    winVec = W(:);
end

function w = local_window(N, typeStr)
    switch lower(typeStr)
        case 'hann'
            n = (0:N-1)';
            if N == 1
                w = 1;
            else
                w = 0.5 - 0.5*cos(2*pi*n/(N-1));
            end
        case 'hamming'
            n = (0:N-1)';
            if N == 1
                w = 1;
            else
                w = 0.54 - 0.46*cos(2*pi*n/(N-1));
            end
        case 'none'
            w = ones(N,1);
        otherwise
            error('Unknown window type: %s', typeStr);
    end
    w = w(:);
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

function out = ternary_str(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
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