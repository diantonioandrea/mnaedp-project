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

	% L2 error (velocity).
	nexttile;

	loglog(sizes, l2Errors, "DisplayName", "Numerical error", ...
		"LineWidth", 2);
	hold on;

	coeffs = polyfit(log(sizes), log(l2Errors), 1);
	loglog(sizes, exp(coeffs(2)) * (sizes .^ coeffs(1)), ...
		"r--", "DisplayName", "Polynomial interpolation", ...
		"LineWidth", 1.5);

	loglog(sizes, sizes, ...
		"g--", "DisplayName", "Linear slope", ...
		"LineWidth", 1.5);
	loglog(sizes, sizes .^ 2, ...
		"b--", "DisplayName", "Quadratic slope", ...
		"LineWidth", 1.5);
	loglog(sizes, sizes .^ 3, ...
		"y--", "DisplayName", "Cubic slope", ...
		"LineWidth", 1.5);

	fprintf("\nVelocity, L2 Error fit: %.2f\n", coeffs(1));
	
	title("L^2 Error: Velocity")

	xlabel("h");
	ylabel("L^2 error")
	
	legend("Location", "northwest");

	hold off;

	% H1 error (velocity).
	nexttile;

	loglog(sizes, h1Errors, "DisplayName", "Numerical error", ...
		"LineWidth", 2);
	hold on;

	coeffs = polyfit(log(sizes), log(h1Errors), 1);
	loglog(sizes, exp(coeffs(2)) * (sizes .^ coeffs(1)), ...
		"r--", "DisplayName", "Polynomial interpolation", ...
		"LineWidth", 1.5);

	loglog(sizes, sizes, ...
		"g--", "DisplayName", "Linear slope", ...
		"LineWidth", 1.5);
	loglog(sizes, sizes .^ 2, ...
		"b--", "DisplayName", "Quadratic slope", ...
		"LineWidth", 1.5);

	fprintf("Velocity, H1 Error fit: %.2f\n", coeffs(1));
	
	title("H^1 Error: Velocity")

	xlabel("h");
	ylabel("H^1 error")
	
	legend("Location", "northwest");
	
	hold off;

	% L2 error (pressure).
	nexttile;

	loglog(sizes, l2ErrorsP, "DisplayName", "Numerical error", ...
		"LineWidth", 2);
	hold on;

	coeffs = polyfit(log(sizes), log(l2ErrorsP), 1);

	loglog(sizes, sizes, ...
		"g--", "DisplayName", "Linear slope", ...
		"LineWidth", 1.5);
	loglog(sizes, sizes .^ 2, ...
		"b--", "DisplayName", "Quadratic slope", ...
		"LineWidth", 1.5);

	fprintf("Pressure, L2 Error fit: %.2f\n", coeffs(1));

	coeffs = polyfit(log(sizes(2:end)), log(l2ErrorsP(2:end)), 1);
	loglog(sizes, exp(coeffs(2)) * (sizes .^ coeffs(1)), ...
		"r--", "DisplayName", "Tail polynomial interpolation", ...
		"LineWidth", 1);

	fprintf("Pressure, L2 Tail error fit: %.2f\n", coeffs(1));
	
	title("L^2 Error: Pressure")

	xlabel("h");
	ylabel("L^2 error")

	legend("Location", "northwest");

	hold off;

	% % Exports plots.
	% exportgraphics(gcf, '../gallery/errorTrend.pdf', ...
	% 		'ContentType', 'vector');
end