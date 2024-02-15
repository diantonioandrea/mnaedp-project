function meshinfo(~)
    fprintf("Stokes CR info.\n\n");

    data = load("quadmeshes.mat");
    input = solver();
	meshes = [data.mesh, data.mesh1, data.mesh2, data.mesh3];

    if nargin > 0
        input.approximate = 1;
		fileID = fopen('../results/infoNB.txt', 'w');
    else
		fileID = fopen('../results/info.txt', 'w');
    end

    for j = 1:4
		fprintf("Solving over mesh: %d\n", j);

        input.mesh = meshes(j);
		output = solver(input);

        fprintf(fileID, "Mesh %d: \n\tSize: %.3f\n\tDOFs: %d\n", ...
            j, output.geometry.size, output.geometry.dofs);
    end
end