% Gradients of basis functions on reference element of [P2 +
% Bubble].

function [gx, gy] = gradBasis(index, x, y)

    switch index
        case 1
            gx = 4 * x + 4 * y - 3;
            gy = 4 * x + 4 * y - 3;
        case 2
            gx = 4 * x - 1;
            gy = 0;
        case 3
            gx = 0;
            gy = 4 * y - 1;
        case 4
            gx = 4 * y;
            gy = 4 * x;
        case 5
            gx = - 4 * y;
            gy =  4 - 8 * y - 4 * x;
        case 6
            gx = 4 - 4 * y - 8 * x;
            gy = -4 * x;

        % Bubble.
        case 7
            gx = 27 * (y - 2 * x * y - y^2);
            gy = 27 * (x - x^2 - 2 * x * y);
    end
end
