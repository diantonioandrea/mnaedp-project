% Errors plotter.
% Andrea Di Antonio, 858798

function errorTrend(~)
	% Arguments.
	narginchk(0, 1);

	tiledlayout(1, 3);
	fprintf("Stokes CR errors.\n\n");

	data = load("quadmeshes.mat");
    input = solver();
	meshes = [data.mesh, data.mesh1, data.mesh2, data.mesh3];

    if nargin > 0
		input.approximate = 1;
    end
	
	l2Errors = zeros(4, 1);
	h1Errors = zeros(4, 1);
	l2ErrorsP = zeros(4, 1);
	
	sizes = zeros(4, 1);

	for j = 1:4
		fprintf("Solving over mesh: %d\n", j);

        input.mesh = meshes(j);
		output = solver(input);

		l2Errors(j) = output.error.l2;
		h1Errors(j) = output.error.h1;
		l2ErrorsP(j) = output.error.l2P;
		
		sizes(j) = output.geometry.size;
	end

	% Colours.
	black = [7 54 66] / 255;
	yellow = [181 137 0] / 255;
	red = [220 50 47] / 255;
	blue = [38 139 210] / 255;
	green = [133 153 0] / 255;

	% Open the file for writing
	if nargin > 0
		fileID = fopen('../results/errorTrendNB.txt', 'w');
	else
		fileID = fopen('../results/errorTrend.txt', 'w');
	end

	% L2 error (velocity).
	nexttile;

	loglog(sizes, l2Errors, "DisplayName", "Numerical error", ...
		"LineWidth", 2, "Color", black);
	hold on;

	coeffs = polyfit(log(sizes), log(l2Errors), 1);
	loglog(sizes, exp(coeffs(2)) * (sizes .^ coeffs(1)), ...
		"--", "DisplayName", "Polynomial interpolation", ...
		"LineWidth", 1, "Color", red);

	% Linear slope.
	shift = l2Errors(1) / sizes(1);
	loglog(sizes, sizes * shift, ...
		"--", "DisplayName", "Linear slope", ...
		"LineWidth", 1, "Color", green);
	
	% Quadratic slope.
	shift = l2Errors(1) / (sizes(1) ^ 2);
	loglog(sizes, sizes .^ 2 * shift, ...
		"--", "DisplayName", "Quadratic slope", ...
		"LineWidth", 1, "Color", blue);

	% Cubic slope.
	shift = l2Errors(1) / (sizes(1) ^ 3);
	loglog(sizes, sizes .^ 3 * shift, ...
		"--", "DisplayName", "Cubic slope", ...
		"LineWidth", 1, "Color", yellow);

	fprintf(fileID, "\nVelocity, L2 Error fit: %.2f\n", coeffs(1));

	title("L^2 Error: Velocity")

	xlabel("h");
	ylabel("L^2 error")

	% H1 error (velocity).
	nexttile;

	loglog(sizes, h1Errors, "DisplayName", "Numerical error", ...
		"LineWidth", 2, "Color", black);
	hold on;

	coeffs = polyfit(log(sizes), log(h1Errors), 1);
	loglog(sizes, exp(coeffs(2)) * (sizes .^ coeffs(1)), ...
		"--", "DisplayName", "Polynomial interpolation", ...
		"LineWidth", 1, "Color", red);

	% Linear slope.
	shift = h1Errors(1) / sizes(1);
	loglog(sizes, sizes * shift, ...
		"--", "DisplayName", "Linear slope", ...
		"LineWidth", 1, "Color", green);
	
	% Quadratic slope.
	shift = h1Errors(1) / (sizes(1) ^ 2);
	loglog(sizes, sizes .^ 2 * shift, ...
		"--", "DisplayName", "Quadratic slope", ...
		"LineWidth", 1, "Color", blue);

	fprintf(fileID, "Velocity, H1 Error fit: %.2f\n", coeffs(1));

	title("H^1 Error: Velocity")

	xlabel("h");
	ylabel("H^1 error")

	% L2 error (pressure).
	nexttile;

	loglog(sizes, l2ErrorsP, "DisplayName", "Numerical error", ...
		"LineWidth", 2, "Color", black);
	hold on;

	coeffs = polyfit(log(sizes), log(l2ErrorsP), 1);

	fprintf(fileID, "Pressure, L2 Error fit: %.2f\n", coeffs(1));

	coeffs = polyfit(log(sizes(2:end)), log(l2ErrorsP(2:end)), 1);
	loglog(sizes, exp(coeffs(2)) * (sizes .^ coeffs(1)), ...
		"--", "DisplayName", "Tail polynomial interpolation", ...
		"LineWidth", 1, "Color", red);

	fprintf(fileID, "Pressure, L2 Tail error fit: %.2f\n", coeffs(1));

	% Linear slope.
	shift = exp(coeffs(2)) * (sizes(1) ^ coeffs(1)) / sizes(1);
	loglog(sizes, sizes * shift, ...
		"--", "DisplayName", "Linear slope", ...
		"LineWidth", 1, "Color", green);
	
	% Quadratic slope.
	shift = exp(coeffs(2)) * (sizes(1) ^ coeffs(1)) / (sizes(1) ^ 2);
	loglog(sizes, sizes .^ 2 * shift, ...
		"--", "DisplayName", "Quadratic slope", ...
		"LineWidth", 1, "Color", blue);

	title("L^2 Error: Pressure")

	xlabel("h");
	ylabel("L^2 error")

	% Close the file
	fclose(fileID);

	% Exports plots.
	plotTitle = '../gallery/errorTrend.pdf';

	if nargin > 0
		plotTitle = '../gallery/errorTrendNB.pdf';
	end

	exportgraphics(gcf, plotTitle, ...
		'ContentType', 'vector', 'BackgroundColor', 'white');

end