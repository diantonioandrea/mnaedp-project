% Loading term at a given coordinate (x,y).

function z = loading(x, y, k)
	x = pi * (x - 0.5);
	y = pi * (y - 0.5);

	switch k
		case 1
			z = (pi^2) * sin(y) * cos(y) * ...
				(1 - 4 * (cos(x))^2) - pi * cos(x);
		case 2
			z = (-pi^2) * sin(x) * cos(x) * ...
				(1 - 4 * (cos(y))^2) + pi * cos(y);
	end
end