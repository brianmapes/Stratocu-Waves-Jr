% Example usage
numPoints = 10;
complexSet1 = randn(numPoints, 1) + 1i * randn(numPoints, 1);
complexSet2 = complexSet1 .* (0.8 + 0.2i) + 0.1 * (randn(numPoints, 1) + 1i * randn(numPoints, 1));

complexScatterPlot(complexSet1, complexSet2);

function complexScatterPlot(complexData1, complexData2)
    % complexData1, complexData2: Arrays of complex numbers.

    [squaredCoherence, meanPhase] = complexScatterAnalysis(complexData1, complexData2);

    figure;
    hold on;

    % Plot complexData1
    scatter(real(complexData1), imag(complexData1), 'filled', 'b'); % Blue markers

    % Plot complexData2
    scatter(real(complexData2), imag(complexData2), 'filled', 'r'); % Red markers

    % Connect corresponding points with thin purple lines
    for i = 1:length(complexData1)
        plot([real(complexData1(i)), real(complexData2(i))], ...
             [imag(complexData1(i)), imag(complexData2(i))], ...
             'Color', [0.5 0 0.5], 'LineWidth', 0.5); % Thin purple lines
    end

    % Add crosshair axes
    xlim_vals = xlim;
    ylim_vals = ylim;
    plot(xlim_vals, [0 0], 'k--', 'LineWidth', 1); % Horizontal axis
    plot([0 0], ylim_vals, 'k--', 'LineWidth', 1); % Vertical axis

    hold off;
    axis equal;
    xlabel('Real');
    ylabel('Imaginary');
    title('Complex Scatter Plot');
    legend('Data 1', 'Data 2');

    % Display squared coherence and mean phase
    text(0.05, 0.95, sprintf('Squared Coherence: %.3f', squaredCoherence), 'Units', 'normalized');
    text(0.05, 0.90, sprintf('Mean Phase: %.3f rad', meanPhase), 'Units', 'normalized');

end

function [squaredCoherence, meanPhase] = complexScatterAnalysis(complexData1, complexData2)
    % complexData1, complexData2: Arrays of complex numbers.

    if length(complexData1) ~= length(complexData2)
        error('Input arrays must have the same length.');
    end

    % Calculate squared coherence
    crossSpectrum = mean(complexData1 .* conj(complexData2));
    powerSpectrum1 = mean(abs(complexData1).^2);
    powerSpectrum2 = mean(abs(complexData2).^2);
    squaredCoherence = abs(crossSpectrum).^2 / (powerSpectrum1 * powerSpectrum2);

    % Calculate mean phase difference
    phaseDifferences = angle(complexData2 ./ complexData1);
    meanPhase = mean(phaseDifferences);
end

