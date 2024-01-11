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
	[xhq, yhq, whq] = quadrature(5);
	Nq = length(xhq);

	% Velocity base size;
	velBase = 7;

	% Basis functions at quadrature nodes on the reference
	% element.
	velPhiQ = zeros(velBase, Nq);

	for i = 1:velBase
    	for q = 1:Nq
        	velPhiQ(i,q) = basis(i, xhq(q), yhq(q));
    	end
	end

	% Pressure basis functions at quadrature nodes on the reference
	% element.
	prePhiQ = zeros(3, Nq);

	for i = 1:3
    	for q = 1:Nq
        	prePhiQ(i,q) = pressureBasis(i, xhq(q), yhq(q));
    	end
	end

	% Gradient of basis functions at quadrature nodes on the
	% reference element.
	gradVelPhiQX = zeros(velBase, Nq);
	gradVelPhiQY = zeros(velBase, Nq);
	
	for i = 1:velBase
    	for q = 1:Nq
        	[gx, gy] = gradBasis(i, xhq(q), yhq(q));
        	gradVelPhiQX(i, q) = gx;
        	gradVelPhiQY(i, q) = gy;
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

		y= [yv(vertex(1)), ...
			yv(vertex(2)), ...
			yv(vertex(3))];

		% Jacobian and derivates.
		jac = [x(2) - x(1), x(3) - x(1); y(2) - y(1), y(3) - y(1)];
		invJac = inv(jac);
		invJacT = invJac';

		% Element's area.
		area = 0.5 * det(jac);

		% Local matrices.
		KE = zeros(2 * velBase, 2 * velBase);
		BE = zeros(3, 2 * velBase);
		LE = zeros(3, 1);
		fE = zeros(2 * velBase, 1);

		% Stiffness matrix quadrature.
		for i = 1:velBase
			for j = 1:i-1
				KE(i, j) = KE(j, i);
			end

			for j = i:velBase
				for q = 1:Nq
					temp = dot((invJacT * [gradVelPhiQX(j, q); gradVelPhiQY(j, q)]), ...
						(invJacT * [gradVelPhiQX(j, q); gradVelPhiQY(j, q)]));

					KE(i, j) = KE(i, j) + nu * temp * whq(q);
				end

				KE(i, j) = 2 * area * KE(i, j);
			end
		end

		KE(velBase + (1:velBase), velBase + (1:velBase)) = KE(1:velBase, 1:velBase);

		% Mixed matrix quadrature.
		for i = 1:(2 * velBase)
			for j = 1:3
				if i < velBase + 0.5
					temp = (invJacT(1, 1:2) * [gradVelPhiQX(i, q); ...
						gradVelPhiQY(i, q)]) * prePhiQ(j, q);
				else
					temp = (invJacT(2, 1:2) * [gradVelPhiQX(i - velBase, q); ...
						gradVelPhiQY(i - velBase, q)]) * prePhiQ(j, q);
				end
			
				BE(j, i) = BE(j, i) + temp * whq(q);
			end

			BE(j, i) = 2 * area * BE(j, i);
		end

		% Loading term quadrature.
		for i = 1:velBase
			for q = 1:Nq
				temp = jac * [xhq(q); yhq(q)] + [x(1); y(1)];

				xq = temp(1);
				yq = temp(2);

				fE(i) = fE(i) + loading(1, xq, yq) * velPhiQ(i, q) * whq(q);
				fE(i + velBase) = fE(i + velBase) + ...
					loading(2, xq, yq) * velPhiQ(i, q) * whq(q);
			end

			fE(i) = 2 * area * fE(i);
			fE(i + velBase) = 2 * area * fE(i + velBase);
		end

		% Multiplier quadrature.
		for i = 1:3
			for q = 1:Nq
				LE(i) = LE(i) + prePhiQ(i, q) * whq(q);
			end

			LE(i) = 2 * area * LE(i);
		end

		% Local to global DOFs.
		velDofsIndexes = [vertex(1), vertex(2), vertex(3), ...
			veNum + edge(1), veNum + edge(2), veNum + edge(3), ...
			veNum + edNum + index];

		doubleVelDofsIndexes = [velDofsIndexes, velDofs + velDofsIndexes];

		preDofsIndexes = [0, 1, 2] * elNum + index;

		% Assembly.
		A(doubleVelDofsIndexes, doubleVelDofsIndexes) = ...
			A(doubleVelDofsIndexes, doubleVelDofsIndexes) + KE;

		B(preDofsIndexes, doubleVelDofsIndexes) = ...
			B(preDofsIndexes, doubleVelDofsIndexes) + BE;
		   
        b(doubleVelDofsIndexes) = b(doubleVelDofsIndexes) + fE;

		L(preDofsIndexes) = L(preDofsIndexes) + LE;
	end

	% BCs.
	if input.approximate == 1
		noBubble = 1:1:(veNum + edNum);
		noBubble = [noBubble, dofs + noBubble];

		A = A(noBubble, noBubble);
		B = B(:, noBubble);
		b = b(noBubble);

		iVertices = setdiff(1:1:veNum, bVertices);
		iEdges = setdiff(1:1:edNum, bEdges);

		intern = [iVertices, veNum + iEdges];
		intern = [intern, veNum + edNum + intern];

		uh = zeros(2 * (veNum + edNum), 1);
	else
		iVertices = setdiff(1:1:veNum, bVertices);
		iEdges = setdiff(1:1:edNum, bEdges);

		intern = [iVertices, veNum + iEdges, veNum + edNum + (1:1:elNum)];
		intern = [intern, velDofs + intern];

		uh = zeros(2 * velDofs, 1);
	end

	A = sparse(A);
	B = sparse(B);

    interns = length(intern);

	Ah = A(intern, intern); clear A;
	Bh = B(:, intern); clear B;
	fh = b(intern); clear b;

	% Global matrix.
	Kh = [Ah, Bh'; Bh, zeros(preDofs)];
	Kh = sparse(Kh);

    % RHS.
    fh = [fh; zeros(preDofs, 1)];
    
	% Multiplier conditions.
	Kh = [Kh, [zeros(interns, 1); L]; ...
		[zeros(1, interns), L', 0]];
	fh = [fh; 0];

	% Linear system solution.
	solh = Kh \ fh;

	% Extracts uh and ph.
	uh(intern) = solh(1:interns);
	uh = reshape(uh, length(uh) / 2, 2);
	ph = solh(interns+1:interns+preDofs);

	% Output.
	output.input = input;
    
    output.matrices.Kh = Kh;
    output.matrices.Ah = Ah;
    output.matrices.Bh = Bh;
    output.matrices.fh = fh;

	output.solution.uh = uh;
	output.solution.ph = ph;
end