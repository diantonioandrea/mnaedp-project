% Pressure basis functions on reference element.

function z = pressurebasis(index, x, y)
	xc = 1/3;
	yc = 1/3;

	switch index
    	case 1
        	z = 1;
    	case 2
        	z = x - xc;
    	case 3
        	z = y - yc;
end