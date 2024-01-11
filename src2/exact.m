% Exact solution at a given (x, y).

function z = exact(index, x, y)
	x = pi * (x - 0.5);
	y = pi * (y - 0.5);

	switch index
		% Velocity, first component.
		case 1
			z = -(cos(x)^2) * cos(y) * sin(y)/2;
		
		% Velocity, second component.
		case 2
			z = (cos(y)^2) * cos(x) * sin(x)/2;

		% Velocity gradient, first component.
		case 3
			z = pi * [cos(x) * sin(x) * cos(y) * sin(y); ...
				-1/2 * (cos(x)^2) * (-sin(y)^2 + cos(y)^2)];
		
		% Velocity gradient, first component.
		case 4
			z = pi * [1/2 * (cos(y)^2) * (-sin(x)^2 + cos(x)^2); ...
				-cos(x) * sin(x) * cos(y) * sin(y)];

		% Pressure.
		case 5
			z = sin(x) - sin(y);
end