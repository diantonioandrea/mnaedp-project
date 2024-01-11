% Errors plotter.
% Andrea Di Antonio, 858798

function errorsPlot(approximate)
	% Arguments.
	narginchk(0, 1);

	if nargin < 1
		approximate = [];
	end

	tiledlayout(1, 3);
	fprintf("Stokes CR errors.\n\n");

	data = load("quadmeshes.mat");
	meshes = [data.mesh, data.mesh1, data.mesh2, data.mesh3];
	
	l2Errors = zeros(4, 1);
	h1Errors = zeros(4, 1);
	l2ErrorsP = zeros(4, 1);
	
	sizes = zeros(4, 1);

	for j = 1:4
		fprintf("Solving over mesh: %d\n", j);

		[~, err, geom] = solver(meshes(j), approximate);

		l2Errors(j) = err.l2;
		h1Errors(j) = err.h1;
		l2ErrorsP(j) = err.l2P;
		
		sizes(j) = geom.diameter;
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

	fprintf("\nVelocity, L2 Error fit: %f\n", coeffs(1));
	
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

	fprintf("Velocity, H1 Error fit: %f\n", coeffs(1));
	
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
	loglog(sizes, exp(coeffs(2)) * (sizes .^ coeffs(1)), ...
		"r--", "DisplayName", "Polynomial interpolation", ...
		"LineWidth", 1.5);

	loglog(sizes, sizes, ...
		"g--", "DisplayName", "Linear slope", ...
		"LineWidth", 1.5);
	loglog(sizes, sizes .^ 2, ...
		"b--", "DisplayName", "Quadratic slope", ...
		"LineWidth", 1.5);

	fprintf("Pressure, L2 Error fit: %f\n", coeffs(1));
	
	title("L^2 Error: Pressure")

	xlabel("h");
	ylabel("L^2 error")

	legend("Location", "northwest");

	hold off;
end