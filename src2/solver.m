% Crouzeix-Raviart element for the Stokes problem.
% Andrea Di Antonio, 858798.

% Adapted from Lorenzo Mascotto's P2-P0 and Mini codes.

function output = solver(input)
	% Solves the Stokes problem given a mesh.

	% Input check.
	narginchk(0, 1);

	% Fixes input.
	if nargin == 0
		data = load("quadmeshes.mat");

		output.mesh = data.mesh;

		output.nu = 1;
		output.approximate = 0;
		output.show = 0;

		return;
	end

	% Mesh informations.
	mesh = input.mesh;

	xv = mesh.xv;
	yv = mesh.yv;

	vertices = mesh.vertices;
	edges = mesh.edges;

	bVertices = mesh.boundary; % Boundary vertices.
	bEdges = mesh.boundedges; % Boundary edges.

	veNum = length(xv); % Vertices' number.
	edNum = size(mesh.endpoints, 1); % Edges' number.
	elNum = size(vertices, 1); % Elements' number.

	% Viscosity constant.
	nu = input.nu;

	% Quadrature informations.
	[xQ, yQ, wQ] = quadrature(5);
	nQ = length(xQ);

	% Velocity base size;
	vBase = 7;

	% Basis functions at quadrature nodes on the reference
	% element.
	vPhiQ = zeros(vBase, nQ);

	for i = 1:vBase
    	for q = 1:nQ
        	vPhiQ(i,q) = basis(i, xQ(q), yQ(q));
    	end
	end

	% Pressure basis functions at quadrature nodes on the reference
	% element.
	pPhiQ = zeros(3, nQ);

	for i = 1:3
    	for q = 1:nQ
        	pPhiQ(i,q) = pressureBasis(i, xQ(q), yQ(q));
    	end
	end

	% Gradient of basis functions at quadrature nodes on the
	% reference element.
	gVPhiQX = zeros(vBase, nQ);
	gVPhiQY = zeros(vBase, nQ);
	
	for i = 1:vBase
    	for q = 1:nQ
        	[gx, gy] = gradBasis(i, xQ(q), yQ(q));

        	gVPhiQX(i, q) = gx;
        	gVPhiQY(i, q) = gy;
    	end
	end

	% DOFs numbers.
	velDofs = veNum + edNum + elNum;
	preDofs = 3 * elNum;

	% Matrices definition.
	
	% Stiffness.
	A = zeros(2 * velDofs);

	% Mixed.
	B = zeros(preDofs, 2 * velDofs);

	% Loading term.
	b = zeros(2 * velDofs, 1);

	% Multiplier.
	L = zeros(preDofs, 1);

	% Assembly.
	for index = 1:elNum
		% Vertices informations.
		vertex = [vertices(index, 1), ...
			vertices(index, 2), ...
			vertices(index, 3)];

		% Edges informations.
		edge = [edges(index, 1), ...
			edges(index, 2), ...
			edges(index, 3)];

		% Coordinates.
		x = [xv(vertex(1)), ...
			xv(vertex(2)), ...
			xv(vertex(3))];

		y = [yv(vertex(1)), ...
			yv(vertex(2)), ...
			yv(vertex(3))];

		% Jacobian and derivates.
		jac = [x(2) - x(1), x(3) - x(1); y(2) - y(1), y(3) - y(1)];
		invJac = inv(jac);
		invJacT = invJac';

		% Element's area.
		area = 0.5 * det(jac);

		% Local matrices.
		sKE = zeros(vBase);
		BE = zeros(3, 2 * vBase);
		LE = zeros(3, 1);
		fE = zeros(2 * vBase, 1);

		% Stiffness matrix quadrature.
		for i = 1:vBase
			for j = 1:i-1
				sKE(i, j) = sKE(j, i);
			end

			for j = i:vBase
				for q = 1:nQ
					temp = dot(invJacT * [gVPhiQX(j, q); gVPhiQY(j, q)], ...
						invJacT * [gVPhiQX(i, q); gVPhiQY(i, q)]);

					sKE(i, j) = sKE(i, j) + nu * temp * wQ(q);
				end

				sKE(i, j) = 2 * area * sKE(i, j);
			end
		end

		% Duplicates scalar stiffness matrix.
		KE = kron(eye(2), sKE);

		% Mixed matrix quadrature.
		for i = 1:(2 * vBase)
			for j = 1:3
				for q = 1:nQ
					if i < (vBase + 0.5)
						temp = (invJacT(1, 1:2) * [gVPhiQX(i, q); ...
							gVPhiQY(i, q)]) * pPhiQ(j, q);
					else
						temp = (invJacT(2, 1:2) * [gVPhiQX(i - vBase, q); ...
							gVPhiQY(i - vBase, q)]) * pPhiQ(j, q);
					end
				
					BE(j, i) = BE(j, i) + temp * wQ(q);
				end

				BE(j, i) = 2 * area * BE(j, i);
			end
		end

		% Loading term quadrature.
		for i = 1:vBase
			for q = 1:nQ
				temp = jac * [xQ(q); yQ(q)] + [x(1); y(1)];

				xq = temp(1);
				yq = temp(2);

				fE(i) = fE(i) + loading(1, xq, yq) * vPhiQ(i, q) * wQ(q);
				fE(i + vBase) = fE(i + vBase) + ...
					loading(2, xq, yq) * vPhiQ(i, q) * wQ(q);
			end

			fE(i) = 2 * area * fE(i);
			fE(i + vBase) = 2 * area * fE(i + vBase);
		end

		% Multiplier quadrature.
		for i = 1:3
			for q = 1:nQ
				LE(i) = LE(i) + pPhiQ(i, q) * wQ(q);
			end

			LE(i) = 2 * area * LE(i);
		end

		% Local to global DOFs.
		vDIndexes = [vertex(1), vertex(2), vertex(3), ...
			veNum + edge(1), veNum + edge(2), veNum + edge(3), ...
			veNum + edNum + index];

		dVDIndexes = [vDIndexes, velDofs + vDIndexes];

		pDIndexes = [0, 1, 2] * elNum + index;

		% Assembly.
		A(dVDIndexes, dVDIndexes) = A(dVDIndexes, dVDIndexes) + KE;
		B(pDIndexes, dVDIndexes) = B(pDIndexes, dVDIndexes) + BE;
        b(dVDIndexes) = b(dVDIndexes) + fE;
		L(pDIndexes) = L(pDIndexes) + LE;
	end

	% BCs.
	if input.approximate == 1
		noBubble = 1:1:(veNum + edNum);
		noBubble = [noBubble, velDofs + noBubble];

		A = A(noBubble, noBubble);
		B = B(:, noBubble);
		b = b(noBubble);

		iVertices = setdiff(1:1:veNum, bVertices);
		iEdges = setdiff(1:1:edNum, bEdges);

		intern = [iVertices, veNum + iEdges];
		
		velDofs = veNum + edNum;
	else
		iVertices = setdiff(1:1:veNum, bVertices);
		iEdges = setdiff(1:1:edNum, bEdges);

		intern = [iVertices, veNum + iEdges, veNum + edNum + (1:1:elNum)];
	end

	intern = [intern, velDofs + intern];

	A = sparse(A);
	B = sparse(B);

	Ah = A(intern, intern); clear A;
	Bh = B(:, intern); clear B;
	fh = b(intern); clear b;

	% Global matrix and RHS.
	Kh = [Ah, Bh'; Bh, zeros(preDofs)];
    fh = [fh; zeros(preDofs, 1)];
    
	% Multiplier conditions.
	interns = length(intern);
	Kh = [Kh, [zeros(interns, 1); L]; ...
		[zeros(1, interns), L', 0]];
	fh = [fh; 0];

	% Linear system solution.
	solh = Kh \ fh;

	% Extracts uh and ph.
	uh = zeros(2 * velDofs, 1);
	uh(intern) = solh(1:interns);
	uh = reshape(uh, length(uh) / 2, 2);
	ph = solh(interns+1:interns+preDofs);

	% Errors evaluation.
	if input.approximate == 1
		vBase = vBase - 1;
	end

	l2Error = 0;
	h1Error = 0;
	l2ErrorP = 0;

	for index = 1:elNum
		% Vertices informations.
		vertex = [vertices(index, 1), ...
			vertices(index, 2), ...
			vertices(index, 3)];

		% Edges informations.
		edge = [edges(index, 1), ...
			edges(index, 2), ...
			edges(index, 3)];

		% Coordinates.
		x = [xv(vertex(1)), ...
			xv(vertex(2)), ...
			xv(vertex(3))];

		y = [yv(vertex(1)), ...
			yv(vertex(2)), ...
			yv(vertex(3))];

		% Jacobian and derivates.
		jac = [x(2) - x(1), x(3) - x(1); y(2) - y(1), y(3) - y(1)];
		invJac = inv(jac);
		invJacT = invJac';

		% Element's area.
		area = 0.5 * det(jac);

		% Local to global DOFs.
		if input.approximate == 1
			vDIndexes = [vertex(1), vertex(2), vertex(3), ...
				veNum + edge(1), veNum + edge(2), veNum + edge(3)];
		else
			vDIndexes = [vertex(1), vertex(2), vertex(3), ...
				veNum + edge(1), veNum + edge(2), veNum + edge(3), ...
				veNum + edNum + index];
		end

		pDIndexes = [0, 1, 2] * elNum + index;
		
		% Solutions at dofs.
		uT1 = uh(vDIndexes, 1);
		uT2 = uh(vDIndexes, 2);
		pT = ph(pDIndexes);

		% Partials.
		l2Err = 0;
		h1Err = 0;
		l2ErrP = 0;

		for q = 1:nQ
			% uh, first and second component, values.
			uh1V = 0;
			uh2V = 0;

			% uh, first and second component, gradient's values.
			uhG1V = [0; 0];
			uhG2V = [0; 0];

			for j = 1:vBase
				uh1V = uh1V + uT1(j) * vPhiQ(j, q);
				uh2V = uh2V + uT2(j) * vPhiQ(j, q);

				uhG1V = uhG1V + uT1(j) * (invJacT * ...
					[gVPhiQX(j, q); gVPhiQY(j, q)]);
				uhG2V = uhG2V + uT2(j) * (invJacT * ...
					[gVPhiQX(j, q); gVPhiQY(j, q)]);
			end

			% ph values.
			phV = pT(1) * pPhiQ(1, q) + ...
				pT(2) * pPhiQ(2, q) + ...
				pT(3) * pPhiQ(3, q);

			% Gauss node on the element.
			node = jac * [xQ(q); yQ(q)] + [x(1); y(1)];

			xq = node(1);
			yq = node(2);

			% Error on node.
			l2Err = l2Err + (exact(1, xq, yq) - uh1V) ^ 2 * wQ(q) + ...
				(exact(2, xq, yq) - uh2V) ^ 2 * wQ(q);

			h1Err = h1Err + norm(exact(3, xq, yq) - uhG1V) ^ 2 * wQ(q) + ...
				norm(exact(4, xq, yq) - uhG2V) ^ 2 * wQ(q);

			l2ErrP = l2ErrP + (exact(5, xq, yq) - phV) ^ 2 * wQ(q);

		end

		l2Err = l2Err * 2 * area;
		h1Err = h1Err * 2 * area;
		l2ErrP = l2ErrP * 2 * area;

		% Total error (squared).
		l2Error = l2Error + l2Err;
		h1Error = h1Error + h1Err;
		l2ErrorP = l2ErrorP + l2ErrP;
	end

	% Mesh size.
	diameter = 0;

	for j = 1:length(mesh.endpoints)
		x1 = mesh.xv(mesh.endpoints(j, 1));
    	y1 = mesh.yv(mesh.endpoints(j, 1));
    	x2 = mesh.xv(mesh.endpoints(j, 2));
    	y2 = mesh.yv(mesh.endpoints(j, 2));

		edge = pdist([x1, y1; x2, y2], "euclidean");

		if edge > diameter
			diameter = edge;
		end
	end

	% Output.

	% Input.
	output.input = input;

	% Matrices.
    output.matrices.Kh = Kh;
    output.matrices.Ah = Ah;
    output.matrices.Bh = Bh;
    output.matrices.fh = fh;

	% Solution.
	output.solution.uh = uh;
	output.solution.ph = ph;

	% Error.
	output.error.l2 = sqrt(l2Error);
	output.error.h1 = sqrt(h1Error);
	output.error.l2P = sqrt(l2ErrorP);

	% Geometry.
	output.geometry.size = diameter;
	output.geometry.dofs = length(fh);

	% Plot, if requested.
	if input.show == 1
		tiledlayout(2, 2);

		% Velocity plot on vertices dofs.
		% Numerical vs analytical.
		nexttile;

		quiver(xv, yv, uh(1:veNum, 1), ...
			uh(1:veNum, 2));

		hold on;
		title("Numerical velocity");
		hold off;
		
		nexttile;

		uex = zeros(veNum, 1);
		uey = zeros(veNum, 1);

		for i = 1:veNum
			uex(i) = exact(1, xv(i), yv(i));
			uey(i) = exact(2, xv(i), yv(i));
		end

		quiver(xv, yv, uex, uey);

		hold on;
		title("Analytical velocity");
		hold off;

		% Pressure plot on centroids.
		nexttile;

		for k = 1:elNum
			hold on;

			indexes = vertices(k, 1:3)';

			pressure = ph(k) * ...
				ones(length(indexes) + 1, 1);
			
			triangle = [xv(indexes), yv(indexes); ...
				xv(indexes(1)), yv(indexes(1))];

			fill3(triangle(:, 1), triangle(:, 2), ...
				pressure, pressure);
		end

		view(2);
		grid on;
		colorbar;
		title("Numerical pressure");
		hold off;

		nexttile;

		for k = 1:elNum
			hold on;

			indexes = vertices(k, 1:3)';
			
			xc = sum(xv(indexes)) / 3;
			yc = sum(yv(indexes)) / 3;

			pressure = exact(5, xc, yc) * ...
				ones(length(indexes) + 1, 1);
			
			triangle = [xv(indexes), yv(indexes); ...
				xv(indexes(1)), yv(indexes(1))];

			fill3(triangle(:, 1), triangle(:, 2), ...
				pressure, pressure);
		end

		view(2);
		grid on;
		colorbar;
		title("Analytical pressure");
		hold off;
	end
end