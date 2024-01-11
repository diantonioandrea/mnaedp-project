% Basis functions on reference element of [P2 + Bubble].

function z = basis(index, x, y)
	switch index
    	case 1
        	z = 2 * x^2 + 4 * x * y - 3 * x + 2 * y^2 - 3 * y + 1;
    	case 2
        	z = x * (2 * x - 1);
    	case 3
        	z = y * (2 * y - 1);

    	case 5
        	z = 4 * x * y;
    	case 6
        	z = - 4 * y * (x + y - 1);
    	case 7
        	z = - 4 * x * (x + y - 1);
		
		% Bubble.
		case 4
        z = 27 * x * y * (1 - x - y);
end
