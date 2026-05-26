clc; clear; close all;
% %%%%%%%%%%%%%%%%% Dataset Selection %%%%%%%%%%%%%%%%%
% dataset:
dataset = 2;

c0 = physconst('lightspeed');

if dataset == 1
    % ==============================================================
    % DATASET 1 YANIK
    % Note: TrueImage size 76x306 ---> ref_optical_image 76x306
    % ==============================================================
    dataPath = '';
    dataName = 'rawdata_yanik';

    bbox = [-70 80 -80 70];
    z0   = 280;

    % Load raw data
    rawData = load([dataPath dataName]);
    rawData = rawData.('rawData3D_simple2D');

    % Parameters
    nFFTtime  = 1024;
    nFFTspace = 1024;
    dx = 200/(size(rawData,3)-1);     % mm
    dy = 198/(size(rawData,2)-1);     % mm

    F0 = 77e9;
    FS = 9121e3;
    K0 = 63.343e12;
    tI = 4.5225e-10;

    % Preprocessing
    rawDataFFT = fft(rawData, nFFTtime);
    k0 = round(K0/FS*(2*z0*1e-3/c0 + tI)*nFFTtime);
    sarData = squeeze(rawDataFFT(k0+1,:,:)); % [100X407]
    %sarData = sarData(:, end:-1:1); 

elseif dataset == 2 || dataset == 3
    % ==============================================================
    % DATASET 2 or 3 (3Dried) 
    % Note: TrueImage size 201x401 ---> ref_optical_image 201x401
    % ==============================================================
    dataPath = '';
    if dataset == 2
        dataName = 'rawdata_gun_3dried';
        bbox = [-150 250 -320 80];
    else
        dataName = 'rawdata_knife_3dried';
        bbox = [-160 250 -210 260];
    end
    % Load raw data
    rawData = load([dataPath dataName]);
    rawData = rawData.('Echo');


    %bbox = [-150 250 -320 80];
    %z0   = 280;
    num_sample = 256;
    Nx = 407;                            % The sampling points in the horizontal direction
    Nz = 200;  
 
    % Parameters
    nFFTtime = num_sample;
    nFFTspace = 512;
    dx = 1;     % mm
    dy = 2;     % mm

    F0 = (77 + 1.8)*1e9;
    %FS = 9121e3;
    Fs = 5*1e6; 
    K0 = 70.295e12; 
    tI = 6.2516e-10;

    % Preprocessing
    Sr = fft(rawData,nFFTtime,3);
    %rawDataFFT = fft(rawData, nFFTtime);
    %k0 = round(K0/FS*(2*z0*1e-3/c0 + tI)*nFFTtime);
    %sarData = squeeze(rawDataFFT(k0+1,:,:));
    %sarData = sarData(:, end:-1:1);
    ID_select = 17;      
    R = c0/2*(ID_select/(K0*(1/Fs)*nFFTtime) - tI);
    z0 = R * 1000;
    sarData = squeeze(Sr(:,:,ID_select)).'; 
 
    for ii = 2:2:Nz
        sarData(ii,:) = fliplr(sarData(ii,:));
    end

elseif dataset == 4 || dataset == 5 || dataset == 6
    % ==============================================================
    % DATASET 4, 5, or 6 ours
    % ==============================================================
    dataPath = '';
    if dataset == 4
        dataName = 'rawdata_knife_ours';
        dx = 1; dy = 1;
        z0   = 185;  
        bbox = [-200 200 -200 200];  
    elseif dataset == 5
        dataName = 'rawdata_plier_ours';
        dx = 1; dy = 2;
        z0   = 210; 
        bbox = [-200 200 -200 200];
    else
        dataName = 'rawdata_screw_driver_ours';
        dx = 1; dy = 2;
        z0   = 230;
        bbox = [-200 200 -200 200];
    end
                  
    % Load raw data
    rawData = load([dataPath dataName]);
    rawData = rawData.('adcDataCube');

    % Parameters
    nFFTtime  = 1024;
    nFFTspace = 1024;
    %dx = 200/(size(rawData,3)-1);      % mm
    %dy = 198/(size(rawData,2)-1);      % mm

    F0 = 77e9;
    FS = 5000e3;
    K0 = 70.295e12;
    tI = 4.5225e-10;

    % Preprocessing
    rawDataFFT = fft(rawData, nFFTtime);
    k0 = round(K0/FS*(2*z0*1e-3/c0 + tI)*nFFTtime);
    sarData = squeeze(rawDataFFT(k0+1,:,:));
    %sarData = sarData(:, end:-1:1);
    for ii = 2:2:size(sarData, 1)
        sarData(ii, :) = fliplr(sarData(ii, :));
    end
else
    error('Unknown dataset selection. Choose dataset = 1, 2, 3, 4, 5, or 6.');
end

% %%%%%%%%%%%%%% Reference SAR Image Generation %%%%%%%%%%%%%%%%%%
% --- Reconstruct SAR image using matched filtering algorithm ---
[trueImageX, trueImageY, trueImage] = ...
    imaging_2DMF(sarData, nFFTspace, z0, dx, dy, bbox, F0);

% Normalize SAR image magnitude to [0,1] for visualization and cpselect
I1_mag = mat2gray(abs(trueImage));   % SAR image (FIXED reference)
[A, B] = size(I1_mag);               % Image dimensions

figure;
imagesc(abs(trueImage));
axis xy;
xlabel('Pixel index (X)');
ylabel('Pixel index (Y)');
title('SAR Image (Pixel Coordinates)');
colormap('jet');

figure;
mesh(trueImageX, trueImageY, abs(trueImage), ...
    'FaceColor','interp','LineStyle','none');
view(2);
xlabel('Horizontal X (mm)');
ylabel('Vertical Y (mm)');
title('SAR Image (Physical Coordinates)');
colormap('jet');
%xlim([-150 250]);
%ylim([-320 80]);
%xlim([-160 250]); 
%xlim([-210 260]);

%% %%%%%%%%%%%%%% Optical Image Preprocessing %%%%%%%%%%%%%%%%%%%%%%%%
%{
% Load the real optical image and prepare it for registration
%I0_raw = im2double(rgb2gray(imread('yanik_real_object2.png'))); % Read & convert to grayscale
%I0_raw = im2double(rgb2gray(imread('knife_target_3Dried_1.jpg')));
%I0_raw = im2double(rgb2gray(imread('knife_target_ours.jpg')));
I0_raw = im2double(rgb2gray(imread('plier_target_ours.jpg')));
I0 = flipud(imresize(I0_raw, size(I1_mag)));                    % Resize to SAR size and flip vertically
% The flipped optical image (I0) will be the MOVING image.

%% INTERACTIVE ALIGNMENT (REPLACES HARDCODED POINTS)
fprintf(['Launching CPSELECT...\n',...
         'Select at least 3 matching points between the Optical (left) and SAR (right) images.\n',...
         'Each N-th point in Optical must correspond to the same physical point in SAR.\n',...
         'When finished, close the window (File → Close CP Select) to continue.\n']);
% Launch cpselect for manual point selection.
% The 'Wait' flag pauses script execution until you close the tool.
[movingPts, fixedPts] = cpselect(I0, I1_mag, 'Wait', true);
%}
%-------------------------------------------03/28
cmap = jet(256);
I1_rgb = ind2rgb(gray2ind(I1_mag, 256), cmap);

%I0_raw = im2double(rgb2gray(imread('plier_target_ours.jpg')));
%I0_raw = im2double(rgb2gray(imread('yanik_real_object2.png')));
%I0_raw = im2double(rgb2gray(imread('knife_target_3Dried.jpg')));
I0_raw = im2double(rgb2gray(imread('gun_target_3Dried_mask.jpg')));
%I0_raw = im2double(rgb2gray(imread('knife_target_ours.png')));
%I0_raw = im2double(rgb2gray(imread('plier_target_ours.png')));
%I0_raw = im2double(rgb2gray(imread('screw_driver_target_ours.png')));

I0 = flipud(imresize(I0_raw, size(I1_mag)));

[movingPts, fixedPts] = cpselect(I0, I1_rgb, 'Wait', true);
%-------------------------------------------

% --- Verify that enough points were selected ---
N_pts = size(movingPts, 1);
if N_pts < 3
    error('At least 3 non-collinear points are required for an Affine transformation. Found %d.', N_pts);
end

%% GEOMETRIC TRANSFORMATION AND WARPING
% --- Estimate the transformation from Optical → SAR space ---
modelType = 'affine';              % Use affine model (rotation, scale, shear, translation)
tform01 = fitgeotrans(movingPts, fixedPts, modelType);

% --- Warp the optical image into SAR pixel grid ---
R1 = imref2d(size(I1_mag));        % Define reference frame = SAR grid
fillWhite = 1.0;                   % Background fill value for outside pixels (white)
I0_warp = imwarp(I0, tform01, ...
    'OutputView', R1, ...
    'InterpolationMethod','bilinear', ...
    'FillValues', fillWhite);

% Save the warped optical image so it can be used later as a reference target.
% This allows pixel-wise comparison or optimization against the SAR image.
ref_optical_image = I0_warp;  % Rename for clarity

save('ref_optical_image.mat', 'ref_optical_image');  % Save in MATLAB format
fprintf('✅ Saved aligned optical image as "ref_optical_image.mat" in the current folder.\n');
%% EVALUATION AND VISUALIZATION
% --- Compute pixel-wise geometric alignment error ---
[mx, my] = transformPointsForward(tform01, movingPts(:,1), movingPts(:,2));
errs    = hypot(mx - fixedPts(:,1), my - fixedPts(:,2));
rmsErr  = sqrt(mean(errs.^2));
fprintf('I0→I1 alignment error (pixels): RMS = %.4f (Target < 2.0)\n', rmsErr);

% --- Display registration results side-by-side ---
figure('Name', 'Registration Verification');
subplot(1,3,1);
imshow(I1_mag, []); title('1. SAR Reference (Fixed)'); axis on;
subplot(1,3,2);
imshow(I0, []); title('2. Optical (Moving, Pre-Warp)'); axis on;
subplot(1,3,3);
imshow(I0_warp, []); title('3. Optical Warped → SAR Grid'); axis on;


%% %%%%%%%%%%%%%% Final Aligned Image Comparison %%%%%%%%%%%%%%%%%%%%%
% --- Compare final aligned images in pixel and physical coordinates ---
img0 = I0_warp;        % Aligned optical image (grayscale)
img1 = trueImage;      % Reference SAR image (complex-valued)

figure('Name', 'Final Aligned Comparison Grid');

% (1) Optical image in pixel coordinates
subplot(2,2,1);
imagesc(img0); axis xy;
xlabel('Pixel index (X)'); ylabel('Pixel index (Y)');
title('Optical Image (Aligned Pixel Coordinates)');
colormap gray;

% (2) Optical image in physical coordinates (mm)
subplot(2,2,2);
imagesc(trueImageX, trueImageY, img0); axis xy;
xlabel('Horizontal X (mm)'); ylabel('Vertical Y (mm)');
title('Optical Image (Physical Coordinates)');
colormap gray;

% (3) SAR image in pixel coordinates
subplot(2,2,3);
imagesc(abs(img1)); axis xy;
xlabel('Pixel index (X)'); ylabel('Pixel index (Y)');
title('SAR Image (Pixel Coordinates)');
colormap('jet');

% (4) SAR image in physical coordinates
subplot(2,2,4);
mesh(trueImageX, trueImageY, abs(img1), ...
    'FaceColor','interp','LineStyle','none');
view(2);
xlabel('Horizontal X (mm)'); ylabel('Vertical Y (mm)');
title('SAR Image (Physical Coordinates)');
colormap('jet');
%xlim([-70, 80]); ylim([-80, 70]);
xlim([-150 250]);
ylim([-320 80]);