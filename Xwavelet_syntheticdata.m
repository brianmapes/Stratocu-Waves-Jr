% Test cross wavelet methods on synthetic data, writes figures as png also

clear all  % All figures and data cleared out 
close all 

%----------- FOLDER/PATH SETTINGS ---------------------------
rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5const_200kmwarpmod2'
rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5const_nowave'
%rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5div_200kmwarpmod2'
%rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5const_400kmwarpmod0'
%rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5div_400kmwarpmod5'

%----------- Scatter plot one particular angle and scale ------
iScatS = 3;
iScatA = 7;

%----------- SPATIAL SCALING & RESIZING --------------------
degrees_per_pixel = 0.04;     % Degrees per pixel (typical for GOES)
km_per_degree = 111.32;       % km per degree
shrinkfactor = 2;             % Image is resized by this factor (2 => half the resolution)
invshrinkfactor = 1 / shrinkfactor;
original_px_km = degrees_per_pixel * km_per_degree;
pixel_size_km = original_px_km * shrinkfactor;

% Seconds between frames (important for speed estimation)
time_resolution = 1800;

%----------- WAVELET PARAMETERS -----------------------------
Angles = 0 : pi/(7*2) : pi;               % 15 Wavelet angles (in radians)
Scales = [2, 4, 8, 16, 32, 64, 128];      % 7 Wavelet scales in pixel units

% Higher res, slower computations
%Angles = 0 : pi/(7*3) : pi;               %22 Wavelet angles (in radians)
%Scales = 2.^( (1:21)/3. )                 %21 scales

NANGLES = numel(Angles);                  % Number of angles
NSCALES = numel(Scales);                  % Number of scales

% For the elongated-envelope wavelet instrad of basic Cauchy, set "true"
CustomWavelet = false;  % This flag indicates whether to use a custom (elliptical) wavelet instead of the default built-in one.
coneAngle = pi/6;      % The variable coneAngle specifies the angular extent of the directional mask and influences the wavelet's sensitivity to orientation.
sigmaX = 0.05;         % The parameter sigmaX defines the decay rate of the wavelet's envelope along the horizontal frequency axis (ωX), affecting its horizontal resolution.
sigmaY = 1.95;         % The parameter sigmaY defines the decay rate of the wavelet's envelope along the vertical frequency axis (ωY), affecting its vertical resolution.
alpha = 0.5;           % The variable alpha is an overall radial decay factor that controls how sharply the wavelet decays in the frequency domain, thereby influencing the scale (frequency) resolution independently of the directional parameters.

%----------- WINDOWING SETTINGS -----------------------------
doWindow = true;              % Flag to apply windowing around image
windowType = 'rectangular';   % 'radial' or 'rectangular'
radius_factor = 0.6;          % Parameter for window function (if used)
decay_rate = 10;              % Controls steepness of window edge

%----------- "ROI" SUBSQUARES-PARTITIONING PARAMETERS ----------------
window_buffer = 0;          % Number of pixels to ignore at each edge
square_size_deg = 5;        % Square size (in degrees) for ROI partitioning


%----------- SMOOTH WAVE-ROSE DISPLAYS, & PEAK DETECTION ---------------------
nAngles_fineFactor  = 4;            % Factor to refine angular resolution in the rose plot
nScales_fineFactor  = 4;            % Factor to refine scale resolution in the rose plot
peakDetectionFactor = 1;            % Threshold factor (mean + factor*std) for peak detection
contourArray        = [95 97 99];   % [Used as either percentiles or absolute values for contouring]
ArrayMode           = 'percentile'; % 'percentile' or 'absolute'

%----------- IMAGE ANNOTATIONS & OUTPUT ---------------------
saverose = true;            % Flag to save the wave–rose image
DisplayValuePower = 1 *10^-4;
DisplayValueCoherence = 1;
DisplayValueSpeed = [-30;30];

 %% GATHER ALL PNG INPUT IMAGE FILENAMES IN DIRECTORY
 if ~exist(rootSepacDir, 'dir')
     error('Folder does not exist: %s', rootSepacDir);
 end
 pngFiles = dir(fullfile(rootSepacDir, '*.png'));
 if isempty(pngFiles)
     fprintf('No PNG images found in %s\n', rootSepacDir);
     return;
 end

numFrames = numel(pngFiles);
fprintf('Found %d frames in dir %s.\n', numFrames);

%% Initializations of empty containers 

prevWaveletSpec = [];    % empty container for swapping in pairwise loop

% Initialize pooling variables for spectrum and xspectrum averages.
% A sum to add up, a count to increment - average is computed from them.
% a "batch" is all the files in the directory, plan accordingly.

batchFrameCount = 0;  % Number of frames processed in loop
batchPairCount = 0;   % Number of pairs processed in loop, 1 less than above
specBatchSum = [];    % Sum for complex spectrum 
powerBatchSum = [];   % Sum for complex spectrum squared amplitude
xspecBatchSum = [];   % Sum for complex cross spectrum 


%% BIG OUTER LOOP OVER FRAMES 
for f_idx = 1:numFrames
    batchFrameCount = batchFrameCount + 1;

%%% GRAB DATA & WINDOW IT
    thisFileName = pngFiles(f_idx).name;
    thisFullPath = strcat(pngFiles(f_idx).folder, '/', pngFiles(f_idx).name);
    fprintf('\nFrame [%d/%d]: %s\n', f_idx, numFrames, thisFullPath);
    
    % READ THE PNG as double precision (range ~[0,1] if 8-bit):
    rawImg = imread(thisFullPath);
    if ndims(rawImg) == 3
        % If it's RGB, convert to grayscale
        rawImg = rgb2gray(rawImg);
    end
    data = double(rawImg) / 255;  % scale to [0..1], or remove /255 if your images are already [0..1]
    data_filt = data;  % Strip all filt/preproc out for tests, data_filt = data
    
    % Resize if needed, turning data_filt into data_pre for the math
    if shrinkfactor ~= 1
        data_pre = imresize(data_filt, invshrinkfactor);
    end
    
    % Apply windowing if enabled.
    if doWindow
        switch lower(windowType)
            case 'radial'
                data_pre = applyRadialWindow(data_pre, radius_factor, decay_rate);
                data_filt = applyRadialWindow(data_filt, radius_factor, decay_rate);
            case 'rectangular'
                data_pre = applyRectangularWindow(data_pre, radius_factor, decay_rate);
                data_filt = applyRectangularWindow(data_filt, radius_factor, decay_rate);
            otherwise
                warning('Unknown window type: %s. No window applied.', windowType);
        end
    end    
    [rowsF, colsF] = size(data_pre);

% Quick look at windowing for sanity check, 0.6 and 10 looks okay
%    if (f_idx == 1)
%        figure; 
%        imshow( circshift(data_filt,[2*max(Scales),2*max(Scales)]));
%        title('Edge artifacts on cared-about scales?');
%    end

%%% COMPUTE 2D COMPLEX TRANSFORM, CALIBRATE TO PROJECTION UNIT (factor 2/S) 
    if CustomWavelet
        waveStruct = barebonesCauchy2D_Elliptical_NoShift(data_pre, Scales, Angles, ...
        coneAngle, sigmaX, sigmaY, alpha);
    else
        waveStruct = cwtft2(data_pre, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    end

    % Empirically, the projection coefficient we want has S/2 factor.
    % Remove that factor so the projection of a unit pure signal is 1. 
    spec_full = squeeze(waveStruct.cfs);  % Dimensions: [Ny, Nx, NSCALES, NANGLES]
    for iS = 1:NSCALES
        spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * (2/Scales(iS));
    end
    
%%% POOLING SPEC AND POWER, first time through just set them ----
    if isempty(specBatchSum)
        specBatchSum = spec_full;
        powerBatchSum= abs(spec_full).^2;
        xspecBatchSum= spec_full*0.;   % right sized container for sums
        xspecBatchSum_nonadv= spec_full*0.;   % right sized container for sums

    else
        specBatchSum = specBatchSum + spec_full;
        powerBatchSum= powerBatchSum+ abs(spec_full).^2; 
    end

%%% POOLING XSPEC IF IT EXISTS ----
    if ~isempty(prevWaveletSpec)
        % Compute the complex cross–wavelet product between previous and current frame.
        crossSpec = prevWaveletSpec .* conj(spec_full);
        xspecBatchSum        = xspecBatchSum + crossSpec;
        batchPairCount= batchFrameCount + 1;
    end

%% PLOTS OF THE RUNNING MEAN RESULT FROM BATCH SO FAR 
    if(f_idx > 1)

% PIV pair 
        velocityField = calculatePIV(prevImage, data, 1);  % 1 is shrinkfactor
        uclip = clip(velocityField(:,:,3),-0.1,0.1);
        vclip = clip(velocityField(:,:,4),-0.1,0.1);
        div = divergence(uclip,vclip);
        skip=3; % display with coarser resolution 
        fun = @(block_struct) mean(block_struct.data, 'all');
        figure; imagesc(blockproc(div,[skip skip],fun),[-1 1]*1e-1); 
        colorbar; hold on; 
        quiver(uclip(1:skip:end,1:skip:end), vclip(1:skip:end,1:skip:end)) 
    
% Roses from batch so far   
        areaxspec = mean(mean(xspecBatchSum, 1, 'omitnan'), 2, 'omitnan')/(batchPairCount);
        areapspec = mean(mean(powerBatchSum, 1, 'omitnan'), 2, 'omitnan')/(batchFrameCount);

% speed from pdif
        cohrose = squeeze( abs(areaxspec)./areapspec );
        pdifrose= squeeze( angle(areaxspec)          );
        speedrose = pdifrose*0;          % right sized container
        for ang = 1:length(Angles)       % at each angle, convert to m/s speed
            speedrose(:,ang) = pdifrose(:,ang)/3.1415 .* ...
                (transpose(Scales)*1000*pixel_size_km) /time_resolution;
        end
    
        figure;
        subplot(221)
        contourf(squeeze(cohrose),30); colorbar; 
        title("coh, frame number ",batchFrameCount)
    
        subplot(222)
        imagesc(squeeze(pdifrose), [-pi pi]); colorbar; 
        title('phase diff')
        h = gca; h.YDir = 'normal'; 
    
        subplot(223)
        imagesc(squeeze(speedrose)); colorbar; 
        title('speed (m/s)')
        h = gca; h.YDir = 'normal';

% Estimate pdifrose_adv from a global windspeed WS and direction WD 
% from speedrose scales 1:5
        speed = mean(speedrose(1:5,:),1); speed360 = [speed -speed(2:end)];

        [WS,WD] = max(speed360);  % Poor mans harmonic fit

% With that estimate, construct an expected pdifrose_adv & subtract it
        pdifrose_adv = pdifrose*0;       % right sized container
        for ang = 1:length(Angles)   % at each angle, convert speed2angle
            advspeed = -WS*sin(Angles(ang) + Angles(WD));
            pdifrose_adv(:,ang) = (advspeed*time_resolution)./ ... 
                (Scales*pixel_size_km*1000) *pi;
        end
        subplot(224)
        imagesc(squeeze(pdifrose - pdifrose_adv), [-pi pi]/5.); colorbar; 
        title('phase dif subtracting advxn')
        h = gca; h.YDir = 'normal';


        % Check winds with a plot 100 that accumulates over pairs
        figure (100); hold on; 
        scatter(Angles, speed); 
        scatter(Angles,-WS*sin(Angles + Angles(WD)))
        title('Advecting wind estimate')


% Keep a scatterplot (111) of one special Scale-Angle bin's coh&speed 
% Scatter what points? 5 degree block averages of course
        npix5deg = 555/pixel_size_km;  % 5 degrees blocks 
        fun = @(block_struct) mean(block_struct.data, 'all');

        specPrevBlocks = blockproc(prevWaveletSpec(:,:,iScatS,iScatA), ...
                                        [npix5deg npix5deg], fun);
        specNowBlocks = blockproc(spec_full(:,:,iScatS,iScatA), ...
                                        [npix5deg npix5deg], fun);
        figure; hold on; 
        complexScatterPlot(specPrevBlocks(:), ... 
                           specNowBlocks (:) );
        figure(111); hold on; 
        complexScatterPlot(specPrevBlocks(:), ... 
                           specNowBlocks (:) );

    end % frame index >1 
    
% Update 'previous' for pairwise loop.
    prevWaveletSpec = spec_full;
    prevImage = data;

end % FRAME LOOP 

fprintf('\nAll done. Single–frame and cross–temporal wavelets complete.\n');

%% COARSE-GRAINING INTO ROI OR SQUARES, ABOUT 5 DEGREES LAT-LON

npix5deg = 555/pixel_size_km;  % 5 degrees blocks 

fun = @(block_struct) mean(block_struct.data, 'all');
new_matrix = blockproc(data_pre, [npix5deg npix5deg], fun);
% figure; imagesc(new_matrix); title('mean data in 5 degree blocks');
%size(new_matrix);  % 9 9, integer, clipping of extra from end

% Ready to use for xspec, spec, etc. except it only workd on 2d arrays!
% Need to loop using blockproc or imresize

[rows, cols] = size(new_matrix); % from above, size from data_pre 
coarse_5deg_spec = zeros(rows, cols, length(Scales), length(Angles));
for S = 1:length(Scales)
    for A = 1:length(Angles)
        coarse_5deg_spec(:,:,S,A) = blockproc(spec_full(:,:,S,A), [npix5deg npix5deg], fun);
    end
end


%% HELPER FUNCTIONS

function data_win = applyRectangularWindow(data_in, radius_factor, decay_rate)
    % Compute global median of the input data
    median_val = median(data_in(:));
    
    % Prepare coordinate system
    [rows, cols] = size(data_in);
    cx = cols/2;
    cy = rows/2;
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Normalized absolute distances from center
    dx = abs(X - cx) / cx;
    dy = abs(Y - cy) / cy;
    
    % R is the "rectangular" distance metric, i.e. the maximum of dx, dy
    R = max(dx, dy);
    
    % Compute window values using a logistic function
    window = 1 ./ (1 + exp(decay_rate * (R - radius_factor)));
    
    % Blend data_in with the median (instead of fading to zero):
    data_win = window .* data_in + (1 - window) .* median_val;
end
%--------------------------------------------------------------------------
function velocityField = calculatePIV(image1, image2, shrinkfactor)
    % Calculate PIV between two images
    [xtable, ytable, utable, vtable, ~, ~, ~] = piv_FFTmulti(...
    image1, image2, ...          % Input images
    64, ...                      % Initial interrogation area (Pass 1)
    32, ...                      % Initial step (Pass 1)
    2, ...                       % Subpixel method (2D Gaussian)
    [], ...                      % mask_inpt (no mask)
    [], ...                      % roi_inpt (full image)
    3, ...                       % passes (3 iterations)
    32, ...                      % int2 (Pass 2 IA)
    16, ...                      % int3 (Pass 3 IA)
    0, ...                       % int4 (unused)
    '*linear', ...               % imdeform method
    0, ...                       % repeat (no repeated correlation)
    1, ...                       % mask_auto (disable autocorrelation)
    1, ...                       % do_linear_correlation (enable)
    0, ...                       % do_correlation_matrices (disable)
    0, ...                       % repeat_last_pass (no)
    0.025 ...                    % delta_diff_min (convergence threshold)
);

    % Downsample the velocity vectors by the shrinkfactor
    downsample_factor = shrinkfactor;
    xtable_ds = xtable(1:downsample_factor:end, 1:downsample_factor:end);
    ytable_ds = ytable(1:downsample_factor:end, 1:downsample_factor:end);
    utable_ds = utable(1:downsample_factor:end, 1:downsample_factor:end);
    vtable_ds = vtable(1:downsample_factor:end, 1:downsample_factor:end);

    % Combine the downsampled velocity vectors into a single matrix
    velocityField = cat(3, xtable_ds, ytable_ds, utable_ds, vtable_ds);
end
%--------------------------------------------------------------------------
function [squaredCoherence, meanPhase] = complexScatterAnalysis(complexData)
    % complexData: A matrix where each column represents a set of complex numbers.

    numSets = size(complexData, 2);
    if numSets < 2
        error('At least two sets of complex numbers are required.');
    end

    % Calculate squared coherence (for the first two columns)
    crossSpectrum = mean(complexData(:, 1) .* conj(complexData(:, 2)));
    powerSpectrum1 = mean(abs(complexData(:, 1)).^2);
    powerSpectrum2 = mean(abs(complexData(:, 2)).^2);
    squaredCoherence = abs(crossSpectrum).^2 / (powerSpectrum1 * powerSpectrum2);

    % Calculate mean phase difference (for the first two columns)
    phaseDifferences = angle(complexData(:, 2) ./ complexData(:, 1));
    meanPhase = mean(phaseDifferences);
end
%----------------------------------------------------------------------