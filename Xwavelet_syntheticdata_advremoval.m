% Test cross wavelet methods on synthetic data, writes figures as png also

clear all  % All figures and data cleared out 
close all 

%----------- FOLDER/PATH SETTINGS ---------------------------
rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5const_200kmwarpmod2';
%rootSepacDir = '/Users/bmapes/Github/stratocu_waves/DATA/syn_uv5const_nowave';
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
Scales = [2, 4, 8, 16, 32]; %5Scales , 64, 128];      % 7 Wavelet scales in pixel units

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


%% COMPUTE 2D COMPLEX TRANSFORM: Normal or Custom  
    waveStruct = cwtft2(data_pre, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec_full = squeeze(waveStruct.cfs);  % Dimensions: [Ny, Nx, NSCALES, NANGLES]
    power_full = abs(spec_full).^2; 
    power_prev = power_full; 

    % Empirically, the projection coefficient we want has S/2 factor.
    % Remove that factor so the projection of a unit pure signal is 1. 
    for iS = 1:NSCALES
        spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * (2/Scales(iS));
    end
    
%%% First time through create containers for BatchSum variables ----
    if isempty(specBatchSum)
        specBatchSum = spec_full;             % initial spectrum   
        powerBatchSum= power_full;            % initial spectrum power
        xspecBatchSum= spec_full*0.;          % right sized container for sums
        xspecBatchSum_nonadv= spec_full*0.;   % right sized container for sums
    end

%% FOR 2nd through end of images and spectra, compute cross info ----
    if ~isempty(prevWaveletSpec)

% Cross Spectrum (complex) 
        crossSpec = prevWaveletSpec .* conj(spec_full);
% Update sums 
        specBatchSum = specBatchSum + spec_full;
        powerBatchSum= powerBatchSum+ abs(spec_full).^2; 
        xspecBatchSum= xspecBatchSum + crossSpec;
        batchPairCount= batchFrameCount + 1;

% Area averages (global here): this pair, or batch to date? 
        areaxspec = mean(mean(crossSpec , 1, 'omitnan'), 2, 'omitnan');

        areapower_prev= mean(mean(abs(prevWaveletSpec).^2, 1), 2);
        areapower_now = mean(mean(abs(spec_full)      .^2, 1), 2);

% Coherence and Speedrose from phase dif (global here) 
        coh2rose= squeeze(abs(areaxspec)./ sqrt(areapower_now.*areapower_prev));
        cohrose = sqrt(coh2rose); 
        pdifrose= squeeze( angle(areaxspec)          );
        speedrose = pdifrose*0;          % right sized container
        for ang = 1:length(Angles)       % at each angle, convert to m/s speed
            speedrose(:,ang) = pdifrose(:,ang)/3.1415 .* ...
                (transpose(Scales)*1000*pixel_size_km) /time_resolution;
        end
    
% Estimate windspeed WS and direction WD from speedrose scales 1:5
        speed = mean(speedrose(1:5,:),1); speed360 = [speed -speed(2:end)];
        [WS,WD] = max(speed360);  % Poor mans harmonic fit! 
        WStrue = 7;              % 7 m/s
        WDtrue = 135. *pi/180.;  % from SE  
        % Check winds with a plot 100 that accumulates over pairs
        figure (100); hold on; 
        scatter(Angles, speed); 
        plot(Angles,-WStrue*sin(Angles + WDtrue    ),'black') % true answer
        plot(Angles,-WS*    sin(Angles + Angles(WD)))
        title('Advecting wind estimate')

% pdifrose_adv from pure WS,WD projection on Angle
        pdifrose_adv = pdifrose*0;   % right sized container
        for ang = 1:length(Angles)   % at each angle, convert speed2angle
%            advspeed = -WS*   sin(Angles(ang) + Angles(WD)); % diagnosed
            advspeed = -WStrue*sin(Angles(ang) + WDtrue);     % give it right #s
            pdifrose_adv(:,ang) = (advspeed*time_resolution)./ ... 
                (Scales*pixel_size_km*1000) *pi;
        end

        diffpdif = squeeze(pdifrose - pdifrose_adv); % simple diff of pdif

 %% Advection propagator, in terms of complex fpec_full. Minus sign crucial!
        propagator = exp(-(pdifrose_adv)*1i);
        spec_full_adv = prevWaveletSpec .* reshape( ...
                        propagator,[1, 1, size(propagator)] );
        spec_full_res = spec_full - spec_full_adv;

% diagnose advective phase change for one scale,angle: sign is right
        figure;
        imagesc( real(prevWaveletSpec(:,:,3,11))  ); hold on 
        contour( real( spec_full     (:,:,3,11)), 'black' );
        contour( real( spec_full_adv (:,:,3,11)), 'red' );
        title('Black = actual next spectrum, red = adv. propagator')

       
% diagnose res for wave's scale,angle :  
        figure(7);
        subplot(221)
        imagesc( real(spec_full    (:,:,iScatS,iScatA+1)-prevWaveletSpec(:,:,iScatS,iScatA+1)));
            title('Actual difference'); colorbar;
% Retrieve the color limits
            colorLimits = clim;

        subplot(222)
        imagesc( real(spec_full_adv(:,:,iScatS,iScatA+1)-prevWaveletSpec(:,:,iScatS,iScatA+1)));
            title('Advective difference'); colorbar;
            clim(colorLimits); % Set color limits to match the first subplot

        subplot(223)
        imagesc( real(spec_full    (:,:,iScatS,iScatA+1)-spec_full_adv  (:,:,iScatS,iScatA+1)) )
            title('Residual (Actual-Advective) difference'); colorbar;
            clim(colorLimits); % Set color limits to match the first subplot

        subplot(224)
        %imagesc( real(spec_full    (:,:,iScatS+1,iScatA)-...
        %              spec_full_res(:,:,iScatS+1,iScatA)) );
        %   title('Residual (Actual-Advective) diff for wave S,A'); colorbar;
        %    clim(colorLimits); % Set color limits to match the first subplot

        
        %% RES_idual spectrum cross spectrum with prevWaveletSpec  
        spec_full_res = spec_full - spec_full_adv;
        crossSpec_res = prevWaveletSpec .* conj(spec_full_res);
% RES Area average (global here) 
        areaxspec_res = mean(mean(crossSpec_res        , 1), 2);
        areapower_res = mean(mean(abs(spec_full_res).^2, 1), 2);
    
% RES Roses 
        coh2rose_res= squeeze(abs(areaxspec_res)./ sqrt(areapower_res.*areapower_prev));
        cohrose_res = sqrt(coh2rose_res); 
        pdifrose_res= squeeze( angle(areaxspec_res) );
        speedrose_res = pdifrose_res*0;  % right sized for resXprevious pdif->speed
        diffspeed     = pdifrose_res*0;  % right sized container, simple difference
        for ang = 1:length(Angles)       % at each angle, convert to m/s speed
            speedrose_res(:,ang) = pdifrose_res(:,ang)/3.1415 .* ...
                (transpose(Scales)*1000*pixel_size_km) /time_resolution;
            diffspeed(:,ang) = diffpdif(:,ang)/3.1415 .* ...
                (transpose(Scales)*1000*pixel_size_km) /time_resolution;
        end
    

%% BLOCK AVGS FOR DIAGNOSTIC SCATTERPLOT 111-113 FOR ONE SELECTED SCALE,ANGLE
% Keep a scatterplot (111-112) of one special Scale-Angle bin's coh&speed 
% Scatter what points? 5 degree block averages, not a zillion pixels
        npix5deg = 555/pixel_size_km;  % 5 degrees blocks 
        fun = @(block_struct) mean(block_struct.data, 'all');

        specPrevBlocks = blockproc(prevWaveletSpec(:,:,iScatS,iScatA), ...
                                        [npix5deg npix5deg], fun);
        specNowBlocks  = blockproc(spec_full      (:,:,iScatS,iScatA), ...
                                        [npix5deg npix5deg], fun);
        xspecBlocks    = blockproc(crossSpec      (:,:,iScatS,iScatA), ...
                                        [npix5deg npix5deg], fun);
        xspecBlocks_res= blockproc(crossSpec_res  (:,:,iScatS,iScatA), ...
                                        [npix5deg npix5deg], fun);

        % A 5 degree map of where on the image the big residual is crossSpec_res
        xamplitude_map = abs( xspecBlocks     ); % total xamplitude
        xamplitude_res = abs( xspecBlocks_res ); % residual xamplitude 

%% DISPLAY DIAGNOSTIC SCATTERPLOTS BY BLOCKS
        figure(111); hold on; 
        complexScatterPlot(specPrevBlocks(:), ... 
                           specNowBlocks (:) );
% residual scatterplot 112
        specNowBlocks_res = blockproc(spec_full_res(:,:,iScatS,iScatA), ...
                                        [npix5deg npix5deg], fun);
        figure; hold on; % this pair scatterplot
        complexScatterPlot(specPrevBlocks(:), ... 
                           specNowBlocks_res (:) );
        figure(112); hold on; % accumulated scatterplot
        complexScatterPlot(specPrevBlocks(:), ... 
                           specNowBlocks_res (:) );


%% Display cross and residual-cross roses including xamplitude block view
% Total cross spectrum 
        figure;
        subplot(221)
        imagesc(squeeze(cohrose),[0.8 1]); colorbar; hold on;
        h = gca; h.YDir = 'normal';
        contour(squeeze(cohrose),'black');

        title("coh frame number ",f_idx)
    
        subplot(222)
        imagesc(squeeze(pdifrose)); colorbar; 
        title('phase diff')
        h = gca; h.YDir = 'normal'; 
    
        subplot(223)
        imagesc(squeeze(speedrose)); colorbar; 
        title('speed (m/s)')
        h = gca; h.YDir = 'normal';

        subplot(224)
        imagesc(xamplitude_map); title("xamplitude in 5deg map"); colorbar
   %     imagesc(squeeze(pdifrose - pdifrose_adv), [-pi pi]/5.); colorbar; 
   %     title('phase dif subtracting advxn')
        h = gca; h.YDir = 'normal';

% Residual xspec after advection removed
        figure;
        subplot(221)
        imagesc(squeeze(cohrose_res),[0 1]); colorbar; hold on;
        h = gca; h.YDir = 'normal';
        contour(squeeze(cohrose_res),'black');
        title("resid coh, frame number ",f_idx)
    
        subplot(222)
        % imagesc(diffpdif); colorbar; % simple diff of phase/speed 
        imagesc(squeeze(pdifrose_res)); colorbar; % What Is This Object? Think more
        title('residual: phase diff')
        h = gca; h.YDir = 'normal'; 
    
        subplot(223)
        % imagesc(squeeze(diffspeed)); colorbar; % simple diff of phase/speed 
        imagesc(squeeze(speedrose_res)); colorbar; % speed of pdif of cross(residual,prev)
        title('residual: speed (m/s)')
        h = gca; h.YDir = 'normal';

        subplot(224) 
        imagesc(xamplitude_res); title("resid xamplitude in 5deg map"); colorbar

% Update all the newly computed sums
        xspecBatchSum_nonadv= xspecBatchSum_nonadv + crossSpec_res;

    end % if ~isempty(prevWaveletSpec)
    
% Update 'previous' spectrum for pairwise loop.
    prevWaveletSpec = spec_full;
    %prevPowerSpec = power_full;
    prevImage = data;

end % FRAME LOOP 

fprintf('\nAll done. Single–frame and cross–temporal wavelets complete.\n');




%% UTIL: COARSE-GRAINING INTO ROI OR SQUARES, ABOUT 5 DEGREES LAT-LON

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
