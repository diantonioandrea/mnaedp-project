% Exact solution at a given (x, y).

function z = exact(x, y, index)
	x = pi * (x - 0.5);
	y = pi * (y - 0.5);

	switch index
		case 1 % Velocity, first component.
			z = -(cos(x)^2) * cos(y) * sin(y)/2;
		case 2 % Velocity, second component.
			z = (cos(y)^2) * cos(x) * sin(x)/2;
		case 3 % Pressure.
			z = sin(x) - sin(y);

		% Velocity gradient.
		case 4 % First component.
			z = pi * [cos(x) * sin(x) * cos(y) * sin(y); ...
				-1/2 * (cos(x)^2) * (-sin(y)^2 + cos(y)^2)];

		case 5 % Second component.
			z = pi * [1/2 * (cos(y)^2) * (-sin(x)^2 + cos(x)^2); ...
				-cos(x) * sin(x) * cos(y) * sin(y)];
end