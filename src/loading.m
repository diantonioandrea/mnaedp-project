% Loading term at a given coordinate (x,y).

function z = loading(x,y,k)
	if k == 1
		z = (pi^2) * sin(pi * (y - 0.5)) * cos(pi * (y - 0.5)) * ...
			(1 - 4 * (cos(pi * (x - 0.5)))^2) - pi * cos(pi * (x - 0.5));
	elseif k == 2
		z = (-pi^2) * sin(pi * (x - 0.5)) * cos(pi * (x - 0.5)) * ...
			(1-4 * (cos(pi * (y - 0.5)))^2) + pi * cos(pi * (y - 0.5));
	end
end