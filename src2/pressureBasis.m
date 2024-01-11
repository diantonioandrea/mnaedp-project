% Pressure basis functions on reference element.

function z = pressureBasis(index, x, y)
	
	switch index
    	case 1
        	z = 1;
    	case 2
        	z = x - 1/3;
    	case 3
        	z = y - 1/3;
	end
end