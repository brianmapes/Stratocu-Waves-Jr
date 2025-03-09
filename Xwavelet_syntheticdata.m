% Test cross wavelet methods on synthetic data, writes figures as png also

clear all  % All figures and data cleared out 
close all 

%----------- FOLDER/PATH SETTINGS ---------------------------
rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5const_200kmwarpmod2'
% rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5const_nowave'

outDir = '/Users/bmapes/Downloads/Results';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

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
Angles = 0 : pi/(7*2) : pi;               % Wavelet angles (in radians)
Scales = [2, 4, 8, 16, 32, 64, 128];      % Wavelet scales in pixel units
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
fprintf('Found %d frames in dir %s.\n', numFrames, outDir);

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

% Quick look at windowing for sanity check 
    if (f_idx == 1)
        figure; imshow(data_filt); title('Windowed data example')
    end

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
    else
        specBatchSum = specBatchSum + spec_full;
        powerBatchSum= powerBatchSum+ abs(spec_full).^2; 
    end

%%% POOLING XSPEC IF IT EXISTS ----
    if ~isempty(prevWaveletSpec)
        % Compute the complex cross–wavelet product between previous and current frame.
        crossSpec_product = prevWaveletSpec .* conj(spec_full);
        xspecBatchSum = xspecBatchSum + crossSpec_product;
        batchPairCount= batchFrameCount + 1;
    end
    
    % Update 'previous' for pairwise loop.
    prevWaveletSpec = spec_full;

%%% PLOTS OF THE RUNNING MEAN RESULT FROM BATCH SO FAR 

    %% FROM POOL MEAN, MAKE ALL-AREA MEAN COH AND PHASE & SPEED ROSES  

    areaxspec = mean(mean(xspecBatchSum, 1, 'omitnan'), 2, 'omitnan')/(batchPairCount);
    areapspec = mean(mean(powerBatchSum, 1, 'omitnan'), 2, 'omitnan')/(batchFrameCount);

    cohrose = abs(areaxspec)./areapspec;
    pdifrose= angle(areaxspec);
    %speed = pdif /3.1415 *Scales *pixel_size_km /time_resolution
    
    figure;
    subplot(121)
    contourf(squeeze(cohrose),30); colorbar; 
    title("coh, frame number ",batchFrameCount)

    subplot(122)
    imagesc(squeeze(pdifrose)); colorbar; 
    title('phase diff')
    h = gca; h.YDir = 'normal'; hold on; 

end % FRAME LOOP 

fprintf('\nAll done. Single–frame and cross–temporal wavelets complete.\n');


%% COARSE-GRAINING INTO ROI OR SQUARES, ABOUT 5 DEGREES LAT-LON

% Initialize a cell array.
%powerBatchSumCells = cell(num_squares_y, num_squares_x);


%% HELPER FUNCTIONS
%==========================================================================

function DisplayAggregatedWaveRose(labelStr, waveRoseMat, Scales, Angles, axHandle, maxVal)
% DISPLAYAGGREGATEDWAVEROSE Displays an aggregated wave-rose plot in a provided axes.
%
%   DISPLAYAGGREGATEDWAVEROSE(labelStr, waveRoseMat, Scales, Angles, axHandle, maxVal)
%
%   Inputs:
%     labelStr    - A string label (e.g., 'Coherence_Square') that will be used 
%                   to set the title and determine the colormap.
%     waveRoseMat - A 2D matrix of size [nScales x nAngles] representing the aggregated
%                   power (or coherence) spectrum.
%     Scales      - A vector of scales.
%     Angles      - A vector of angles (in radians).
%     axHandle    - A handle to an existing axes in which to display the wave-rose.
%     maxVal      - (Optional) A scalar maximum value to set the color axis to [0 maxVal].
%
%   Example:
%     % Suppose axInset is an inset axes handle and maxVal is computed globally.
%     DisplayAggregatedWaveRose('Coherence_Square', coherenceMatrix, Scales, Angles, axInset, maxVal);

    % Check that axHandle is provided.
    if nargin < 5 || isempty(axHandle)
        error('An axes handle (axHandle) must be provided.');
    end
    
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        minVal = maxVal(1);
        maxVal = maxVal(2);
    end

    % Set the current axes to the provided handle and hold it.
    axes(axHandle);
    hold(axHandle, 'on');

    % Determine dimensions.
    nScales = length(Scales);
    nAngles = length(Angles);
    
    % Use the provided matrix as the "innerpower" for display.
    innerpower = waveRoseMat;
    
    % Refine grid for smoother plotting.
    nAngles_fineFactor = 4;
    nScales_fineFactor = 4;
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor * nAngles);
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor * nScales);
    
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);
    
    % Interpolate the innerpower onto the fine grid.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';

    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        % For each row (corresponding to a fine-scale value), multiply by that scale.
        for iRow = 1:size(innerpower_fine,1)
            innerpower_fine(iRow, :) = (8900*(innerpower_fine(iRow, :)/(2*pi)) .* Scales_fine_linear(iRow)* pi/sqrt(2))/ 1800;
        end
    end
    
    % Convert polar coordinates to Cartesian (using a log-scale for R).
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);
    
    % Plot the positive half of the rose.
    pcolor(axHandle, X_pos_fine, Y_pos_fine, innerpower_fine);
    shading(axHandle, 'interp');
    
    % Choose colormap based on label.
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        colormap(axHandle, 'jet');
    elseif startsWith(labelStr, 'Coherence', 'IgnoreCase', true)
        colormap(axHandle, 'turbo');
    else
        colormap(axHandle, 'parula');
    end
    axis(axHandle, 'equal', 'tight', 'off');
    
    % Optionally set color axis limits if maxVal is provided.
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
            clim(axHandle, [minVal, maxVal]);
    elseif nargin >= 6 && ~isempty(maxVal)
            clim(axHandle, [0, maxVal]);
    end
    
    % Optionally add radial grid lines.
    max_scale_log = log10(max(Scales));
    for i = 1:length(Scales)
        level_log = log10(Scales(i));
        theta_ring = linspace(0, 2*pi, 180);
        [x_ring, y_ring] = pol2cart(theta_ring, level_log);
        plot(axHandle, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        %text(axHandle, level_log, 0, sprintf('%.1f', Scales(i)), ...
        %     'Color','k','FontSize',4,'HorizontalAlignment','left');
    end

    
    % Add angular lines.
    angle_ticks = linspace(0, 2*pi, 7);
    max_r = max_scale_log * 1.1;
    for i = 1:length(angle_ticks)
        [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
        line(axHandle, [0 x_label], [0 y_label], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    end
end

function produceAggregatedWaveRose(labelStr, waveRoseMat, Scales, Angles, outDir, outName, saverose, maxVal)
    if nargin < 7
        saverose = true;
    end
    if nargin < 8
        maxVal = [];  % if not provided, use auto-scale
    end

    if size(maxVal,1)==2
        minVal = maxVal(1);
        maxVal = maxVal(2);
    end

    % Dimensions
    nScales = length(Scales);
    nAngles = length(Angles);
    
    % Use waveRoseMat as the "innerpower" (power spectrum display)
    innerpower = waveRoseMat;
    
    % Refine grid for smoother plotting.
    nAngles_fineFactor  = 4;
    nScales_fineFactor  = 4;
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor * nAngles);
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor * nScales);
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);
    
    % Interpolate the innerpower onto the fine grid.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';

    % --- Modification for Speed: multiply each row by its corresponding scale ---
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        % For each row (corresponding to a fine-scale value), multiply by that scale.
        for iRow = 1:size(innerpower_fine,1)
            innerpower_fine(iRow, :) = (8900*(innerpower_fine(iRow, :)/(2*pi)) .* Scales_fine_linear(iRow)* pi/sqrt(2))/ 1800;
        end
    end
    
    % Convert to Cartesian coordinates for plotting.
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);
    
    % Create figure.
    figRose = figure('visible','off');
    set(figRose, 'Position', [100 100 600 600]);
    ax1 = axes('Position',[0.1 0.1 0.75 0.75]);
    hold(ax1, 'on');
    pcolor(ax1, X_pos_fine, Y_pos_fine, innerpower_fine);
    shading(ax1, 'interp');
    if startsWith(labelStr, 'Coherence', 'IgnoreCase', true)
        colormap(ax1, 'turbo');
    elseif startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        colormap(ax1,'jet')
    else
        colormap(ax1, 'parula');
    end
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        clim(ax1, [minVal, maxVal]);
    elseif ~isempty(maxVal)  
        clim(ax1, [0, maxVal]);
    else
        clim(ax1, [0, max(innerpower_fine(:))]);
    end
    axis(ax1, 'equal', 'tight', 'off');
    
    % Duplicate plot for negative angles.
    ax2 = axes('Position', ax1.Position, 'Color','none', 'HitTest','off');
    hold(ax2, 'on');
    pcolor(ax2, X_neg_fine, Y_neg_fine, innerpower_fine);
    shading(ax2, 'interp');
    if ~isempty(maxVal)
        clim(ax2, [0, maxVal]);
    else
        clim(ax2, [0, max(innerpower_fine(:))]);
    end
    axis(ax2, 'equal', 'tight', 'off');
    
    uistack(ax1, 'top');
    linkprop([ax1 ax2], {'XLim','YLim','Position','CameraPosition','CameraUpVector'});
    
    % Add radial grid lines (using logarithmic scale for the scales).
    max_scale_log = log10(max(Scales));
    for i = 1:length(Scales)
        level_log = log10(Scales(i));
        theta_ring = linspace(0, 2*pi, 180);
        [x_ring, y_ring] = pol2cart(theta_ring, level_log);
        plot(ax1, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        plot(ax2, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        val_linear = 10^(level_log);
        text(ax1, level_log, 0, sprintf('%.2f', val_linear), 'Color','k','FontSize',4,'HorizontalAlignment','left');
    end
    
    % Add angular lines.
    angle_ticks = linspace(0, 2*pi, 7);
    max_r = max_scale_log * 1.1;
    for i = 1:length(angle_ticks)
        [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
        line(ax1, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
        line(ax2, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
    end
    
    title(ax1, sprintf('%s Wave–Rose', labelStr), 'FontSize',12, 'FontWeight','bold');
    
    c = colorbar(ax1, 'Location','eastoutside');
    c.Label.String = 'Wavelet Power';
    c.Label.FontWeight = 'bold';
    
    ax1_pos = ax1.Position;
    ax1_pos(3) = ax1_pos(3) * 0.85;
    ax1.Position = ax1_pos;
    ax2.Position = ax1_pos;
    
    if saverose
        roseFileName = fullfile(outDir, sprintf('%s.png', outName));
        fprintf('Saving wave–rose plot to: %s\n', roseFileName);
        exportgraphics(figRose, roseFileName, 'Resolution',300);
    end
    close(figRose);
end
%--------------------------------------------------------------------------

function results = extractWaveletFeatures(dataType, spec_full, data_background, squares, ...
    Scales, Angles, frameDateStr, ...
    nAngles_fineFactor, nScales_fineFactor, peakDetectionFactor, contourOption, contourArray)
% EXTRACTWAVELETFEATURES
%  Extracts wavelet features such as the number of peaks, scale/angle pairs,
%  and contour positions for each frame.

    %% 1) Dimensions & Basic Summaries
    [Ny_sh, Nx_sh, nScales, nAngles] = size(spec_full);
    [Ny_orig, Nx_orig] = size(data_background);

    power = abs(spec_full).^2;  % wavelet power
    innerpower = squeeze(mean(mean(power, 1, 'omitnan'), 2, 'omitnan'));

    %% 2) Wave-Rose Visualization with Logarithmic Radial Scale
    % 2.1) Interpolate power to a finer grid
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor*nAngles);
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor*nScales);

    % Create meshgrids for the original and fine grids:
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);

    % Transform the radial coordinate to a logarithmic scale.
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);

    % Interpolate the innerpower on the original (linear) grid and then use
    % the logarithmic radial coordinates for plotting.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';

    % Compute cartesian coordinates for the fine grid using the log-transformed radius:
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);

    % 2.2) Peak detection based on threshold (for labeling)
    threshold_orig = mean(innerpower(:)) + peakDetectionFactor * std(innerpower(:));
    bwMask_orig = (innerpower >= threshold_orig);
    CC = bwconncomp(bwMask_orig, 4);
    numPeaks = CC.NumObjects;

    % Draw contour lines on the coarse (original) polar grid.
    [X_orig, Y_orig] = pol2cart(Theta_orig, R_orig_log);

    % Label each connected region (peak) on the polar plot.
    peakRegions = struct();
    for pk = 1:numPeaks
        [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
        meanScale = mean(Scales(scaleIndices), 'omitnan');
        meanAngle = mean(Angles(angleIndices), 'omitnan');
        peakRegions(pk).ScaleIndices = scaleIndices;
        peakRegions(pk).AngleIndices = angleIndices;
        peakRegions(pk).MeanScale = meanScale;
        peakRegions(pk).MeanAngle = meanAngle;
    end

    %% 3) Region Summaries & Overlays (Final Annotated Image)
    scaleFactorX = Nx_orig / Nx_sh;
    scaleFactorY = Ny_orig / Ny_sh;

    for pk = 1:numPeaks
        waveSum = zeros(Ny_sh, Nx_sh);
        wavePower = zeros(Ny_sh, Nx_sh);
        currentRegion = peakRegions(pk);

        for jj = 1:numel(currentRegion.ScaleIndices)
            s_idx = currentRegion.ScaleIndices(jj);
            a_idx = currentRegion.AngleIndices(jj);
            coeff = spec_full(:,:,s_idx,a_idx);
            waveSum = waveSum + real(coeff);
            wavePower = wavePower + abs(coeff).^2;
        end

        waveSum_up = imresize(waveSum, [Ny_orig, Nx_orig]);
        wavePower_up = imresize(wavePower, [Ny_orig, Nx_orig]);

        % ----- NEW CONTOUR LEVEL SYSTEM -----
        % Choose contour levels based on either absolute values or percentiles.
        switch lower(contourOption)
            case 'absolute'
                % Use the provided absolute values (assumed positive) for contours.
                contourLevels = contourArray;
            case 'percentile'
                % Compute the given percentiles on the absolute values of waveSum_up.
                contourLevels = prctile(abs(waveSum_up(:)), contourArray);
            otherwise
                error('Unknown contour option. Choose either "absolute" or "percentile".');
        end

        % Store contour positions
        contourPositions = struct();
        contourPositions.Positive = contourc(waveSum_up, contourLevels);
        contourPositions.Negative = contourc(waveSum_up, -contourLevels);

        % Store results for the current peak
        peakRegions(pk).ContourPositions = contourPositions;
    end

    % Store results for the current frame
    
    results.NumPeaks = numPeaks;
    results.PeakRegions = peakRegions;
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

function produceAnnotatedImages(dataType, spec_full, data_background, squares, ...
    Scales, Angles, outDir, frameDateStr, ...
    saverose, nAngles_fineFactor, nScales_fineFactor, ...
    peakDetectionFactor,contourOption, contourArray)
% PRODUCEANNOTATEDIMAGES
%  Creates a wave-rose plot (in polar form, duplicated for ± angles),
%  detects peaks based on a threshold, and creates overlay images for each
%  peak region with wavelet real-part and power contours.


    %% 1) Dimensions & Basic Summaries
    [Ny_sh, Nx_sh, nScales, nAngles] = size(spec_full);
    [Ny_orig, Nx_orig] = size(data_background);

    power = abs(spec_full).^2;  % wavelet power
    innerpower = squeeze(mean(mean(power, 1, 'omitnan'), 2, 'omitnan'));

    %% 2) Wave-Rose Visualization with Logarithmic Radial Scale
    % 2.1) Interpolate power to a finer grid
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor*nAngles);
    % Use a logarithmic spacing for scales:
    % First, create a linear grid and then take the logarithm.
    % Alternatively, you could directly use logspace. Here we demonstrate by
    % applying the logarithm to the linear grid.
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor*nScales);
    
    % Create meshgrids for the original and fine grids:
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);
    
    % Transform the radial coordinate to a logarithmic scale.
    % (Assumes Scales > 0.)
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);
    
    % Interpolate the innerpower on the original (linear) grid and then use
    % the logarithmic radial coordinates for plotting.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';
    
    % Compute cartesian coordinates for the fine grid using the log-transformed radius:
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);
    
    figRose = figure('visible','off');
    ax1 = axes('Position',[0.1 0.1 0.75 0.75]);
    hold(ax1, 'on');
    pcolor(ax1, X_pos_fine, Y_pos_fine, innerpower_fine);
    shading(ax1, 'interp');
    colormap(ax1, 'parula');
    axis(ax1, 'equal', 'tight', 'off');
    
    ax2 = axes('Position', ax1.Position, 'Color','none', 'HitTest','off');
    hold(ax2, 'on');
    pcolor(ax2, X_neg_fine, Y_neg_fine, innerpower_fine);
    shading(ax2, 'interp');
    axis(ax2, 'equal', 'tight', 'off');
    
    uistack(ax1, 'top');
    linkprop([ax1 ax2], {'XLim','YLim','Position','CameraPosition','CameraUpVector'});
    
    % 2.2) Peak detection based on threshold (for labeling)
    threshold_orig = mean(innerpower(:)) + peakDetectionFactor * std(innerpower(:));
    bwMask_orig = (innerpower >= threshold_orig);
    CC = bwconncomp(bwMask_orig, 4);
    numPeaks = CC.NumObjects;
    
    % Draw contour lines on the coarse (original) polar grid.
    [X_orig, Y_orig] = pol2cart(Theta_orig, R_orig_log);
    contour(ax1, X_orig, Y_orig, innerpower, [threshold_orig threshold_orig], 'r-', 'LineWidth',2);
    contour(ax2, X_orig, Y_orig, innerpower, [threshold_orig threshold_orig], 'r-', 'LineWidth',2);
    
    % Label each connected region (peak) on the polar plot.
    for pk = 1:numPeaks
        [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
        meanScale = mean(Scales(scaleIndices), 'omitnan');
        meanAngle = mean(Angles(angleIndices), 'omitnan');
        [x_peak, y_peak] = pol2cart(meanAngle, log10(meanScale));
        text(ax1, x_peak, y_peak, sprintf('%d', pk), ...
            'Color','k','FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    %% 2.3) Annotations and Grid Overlays (using logarithmic radial axis)
    % Radial circles: plot circles at each original scale, but use log10(scales)
    for i = 1:length(Scales)
        theta_ring = linspace(0, 2*pi, 100);
        [x_ring, y_ring] = pol2cart(theta_ring, log10(Scales(i)));
        plot(ax1, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        plot(ax2, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        % Label the circle with the original scale value.
        text(ax1, log10(Scales(i))*1.05, 0, sprintf('%.1f', Scales(i)), ...
            'HorizontalAlignment','left','FontSize',4);
    end
    
    % Angular lines
    angle_ticks = linspace(0, 2*pi, 13);
    angle_labels = {'0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6','\pi',...
                    '7\pi/6','4\pi/3','3\pi/2','5\pi/3','11\pi/6','2\pi'};
    max_r = log10(max(Scales)) * 1.1;
    
    for i = 1:length(angle_ticks)
        [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
        line(ax1, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
        line(ax2, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
        if angle_ticks(i) <= pi
            text(ax1, x_label*1.05, y_label*1.05, angle_labels{i}, ...
                'HorizontalAlignment','center','FontSize',8);
        else
            text(ax2, x_label*1.05, y_label*1.05, angle_labels{i}, ...
                'HorizontalAlignment','center','FontSize',8);
        end
    end
    
    % Colorbar for the polar plot
    c = colorbar(ax1, 'Location','eastoutside');
    c.Label.String = 'Wavelet Power';
    c.Label.FontWeight = 'bold';
    ax1_pos = ax1.Position;
    ax1_pos(3) = ax1_pos(3) * 0.85;
    ax1.Position = ax1_pos;
    ax2.Position = ax1_pos;
    
    title(ax1, sprintf('Polar Wave-Rose: %d Significant Regions', numPeaks), ...
        'FontSize',12, 'FontWeight','bold');
    
    if saverose
        roseName = fullfile(outDir, sprintf('WaveRose_%s.png', frameDateStr));
        exportgraphics(figRose, roseName, 'Resolution',300);
    end
    close(figRose);
    
    %% 3) Region Summaries & Overlays (Final Annotated Image)
    scaleFactorX = Nx_orig / Nx_sh;
    scaleFactorY = Ny_orig / Ny_sh;
    
    peakRegions = cell(numPeaks,1);
    for pk = 1:numPeaks
        [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
        scales     = Scales(scaleIndices);
        angles_deg = rad2deg(Angles(angleIndices));
    
        scale_str = join(split(num2str(scales,'%.1f ')), '/');
        angle_str = join(split(num2str(angles_deg,'%.0f ')), '/');
    
        peakRegions{pk} = struct('ScaleIndices',scaleIndices,...
                                 'AngleIndices',angleIndices,...
                                 'ScaleStr',scale_str{1},...
                                 'AngleStr',angle_str{1});
    end
    
    for pk = 1:numPeaks
        waveSum = zeros(Ny_sh, Nx_sh);
        wavePower = zeros(Ny_sh, Nx_sh);
        currentRegion = peakRegions{pk};
    
        for jj = 1:numel(currentRegion.ScaleIndices)
            s_idx = currentRegion.ScaleIndices(jj);
            a_idx = currentRegion.AngleIndices(jj);
            coeff = spec_full(:,:,s_idx,a_idx);
            waveSum = waveSum + real(coeff);
            wavePower = wavePower + abs(coeff).^2;
        end
    
        waveSum_up   = imresize(waveSum,   [Ny_orig, Nx_orig]);
        wavePower_up = imresize(wavePower, [Ny_orig, Nx_orig]);
    
        fig = figure('visible','off');
        switch upper(dataType)
            case 'IR'
                imagesc(data_background, [0 1])
            case 'VIS'
                image(data_background);
            otherwise
                error('Unknown dataType.');
        end
        colormap(gray);
        axis image off;
        hold on;
    
        % ----- NEW CONTOUR LEVEL SYSTEM -----
        % Choose contour levels based on either absolute values or percentiles.
        switch lower(contourOption)
            case 'absolute'
                % Use the provided absolute values (assumed positive) for contours.
                % Draw red contours at the positive levels and blue contours at the corresponding negative levels.
                contourLevels = contourArray;
            case 'percentile'
                % Compute the given percentiles on the absolute values of waveSum_up.
                contourLevels = prctile(abs(waveSum_up(:)), contourArray);
            otherwise
                error('Unknown contour option. Choose either "absolute" or "percentile".');
        end
    
        % Plot contours:
        % For positive values:
        contour(waveSum_up, contourLevels, 'LineColor','red', 'LineWidth',0.5);
        % For negative values (mirror the levels):
        contour(waveSum_up, -contourLevels, 'LineColor','blue', 'LineWidth',0.5);
    
        % ----- Draw ROI squares in final image -----
        for sq = 1:numel(squares)
            xPos_orig = squares(sq).x_range(1) * scaleFactorX;
            yPos_orig = squares(sq).y_range(1) * scaleFactorY;
            w_orig    = length(squares(sq).x_range) * scaleFactorX;
            h_orig    = length(squares(sq).y_range) * scaleFactorY;
            rectangle('Position',[xPos_orig, yPos_orig, w_orig, h_orig],...
                'EdgeColor','k','LineWidth',1);
        end
    
        titleText = {sprintf('Instrument X - %s', frameDateStr), ...
                     sprintf('Peak %d/%d - Scales: %s', pk, numPeaks, peakRegions{pk}.ScaleStr), ...
                     sprintf('Angles: %s°', peakRegions{pk}.AngleStr)};
        title(titleText, 'Color','k','FontWeight','bold','FontSize',10,'Interpreter','none');
    
        outName = fullfile(outDir, sprintf('Frame_%s_Region%02d.png', frameDateStr, pk));
        saveas(fig, outName);
        close(fig);
    end
end
%--------------------------------------------------------------------------
