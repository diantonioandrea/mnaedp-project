% Crouzeix-Raviart element for the Stokes problem.
% Andrea Di Antonio, 858798.

% Adapted from Lorenzo Mascotto's P2-P0 and Mini codes.

function [solution, error, geometry] = solver(mesh, approximate, show)
	% Solves the Stokes problem given a mesh.
	
	% Arguments.
	narginchk(1, 3);

	if nargin < 2
		approximate = [];
	end

	if nargin < 3
		show = [];
	end

	% Extracts geometric information from the given mesh.
	xv = mesh.xv; % Vertices' x.
	yv = mesh.yv; % Vertices' y.
	vertices = mesh.vertices;
	edges = mesh.edges;
	endpoints = mesh.endpoints;
	boundary = mesh.boundary;
	boundedges = mesh.boundedges;

	% Viscosity constant.
	nu = 1;

	% Quadrature nodes.
	[xhq, yhq, whq] = quadrature(5);
	Nq = length(xhq);

	% Constants.
	verticesNumber = length(xv);
	edgesNumber = size(endpoints, 1);
	elementsNumber = size(vertices, 1);

	% Basis functions at quadrature nodes on the reference
	% element.
	phihq = zeros(7, Nq);

	for i = 1:7
    	for q = 1:Nq
        	phihq(i,q) = basis(i, xhq(q), yhq(q));
    	end
	end

	% Pressure basis functions at quadrature nodes on the reference
	% element.
	pphihq = zeros(3, Nq);

	for i = 1:3
    	for q = 1:Nq
        	pphihq(i,q) = pressurebasis(i, xhq(q), yhq(q));
    	end
	end

	% Gradient of basis functions at quadrature nodes on the
	% reference element.
	gphihqx = zeros(7, Nq);
	gphihqy = zeros(7, Nq);
	
	for i = 1:7
    	for q = 1:Nq
        	[gx,gy] = gradbasis(i, xhq(q), yhq(q));
        	gphihqx(i, q) = gx;
        	gphihqy(i, q) = gy;
    	end
	end

	% Matrices.

	% Stiffness matrix.
	A = zeros(2 * (verticesNumber + edgesNumber + elementsNumber));
	
	% Mixed matrix.
	B = zeros(3 * elementsNumber, ...
		2 * (verticesNumber + edgesNumber + elementsNumber));

	% Loading part.
	b = zeros(2 * (verticesNumber + edgesNumber + ... 
		elementsNumber), 1);
	
	% Zero average constraint.
	areaVector = zeros(3 * elementsNumber, 1);

	% Assembly.
	for index = 1:elementsNumber
		% Vertices information.
		vertex = [vertices(index, 1), vertices(index, 2), ...
			vertices(index, 3)];
		
		X = [xv(vertex(1)), xv(vertex(2)), xv(vertex(3))];
		Y = [yv(vertex(1)), yv(vertex(2)), yv(vertex(3))];

		% Jacobian.
		jacob = [X(2) - X(1), X(3) - X(1); Y(2) - Y(1), Y(3) - Y(1)];

		% Inverse Jacobian.
		invJacob = inv(jacob);

		% Transpose inverse Jacobian.
		invJacobT = invJacob';

		% Area.
		area = 0.5 * det(jacob);	

		% Local matrices.
		sKE = zeros(7, 7); % Single KE.
		BE = zeros(3, 14);
		
		% Stiffness matrix quadrature.
		for i = 1:7
			for j = 1:i-1
				sKE(i,j) = sKE(j,i);
			end

        	for j = i:7
				for q = 1:Nq
                	tmp = dot((invJacobT * [gphihqx(j, q); gphihqy(j, q)]),...
                    	(invJacobT * [gphihqx(i, q); gphihqy(i, q)]));

                	sKE(i,j) = sKE(i,j) + nu * tmp * whq(q);
				end

            	sKE(i, j) = 2 * area * sKE(i, j);
        	end
		end

		% Duplicates scalar stiffness matrix.
		KE = kron(eye(2), sKE);

		% Mixed matrix quadrature.
		for i = 1:14
			for j = 1:3
				for q = 1:Nq
					if i < 7.5
                    		tmp = (invJacobT(1, 1:2) * ...
								[gphihqx(i, q); gphihqy(i, q)]) * ...
                        		pphihq(j, q);
            		elseif i > 7.5
                    		tmp = (invJacobT(2, 1:2) * ...
								[gphihqx(i - 7, q); gphihqy(i - 7, q)])*...
                        		pphihq(j, q);
					end
                	
					BE(j,i) = BE(j,i) + tmp * whq(q);
				end

            	BE(j,i) = 2 * area * BE(j,i);
			end
		end

		% Edges information.
		edge = [edges(index, 1), edges(index, 2), edges(index, 3)];

		% 'Local to global' arrays of dofs.
		dofg = [vertex(1), vertex(2), vertex(3), ...
			verticesNumber + index, ...
			verticesNumber + elementsNumber + edge(1), ...
			verticesNumber + elementsNumber + edge(2), ...
			verticesNumber + elementsNumber + edge(3)];

		dofgg = [dofg, dofg + verticesNumber + edgesNumber + ...
			elementsNumber];

		dofp = [0, elementsNumber, 2 * elementsNumber] + index;
		
		% Loading term quadrature.
    	fE = zeros(14, 1);
	
		for i = 1:7
			for q = 1:Nq
            	tmp = jacob * [xhq(q); yhq(q)] + [X(1); Y(1)];

            	xq = tmp(1);
            	yq = tmp(2);

            	fE(i) = fE(i) + loading(xq, yq, 1) * ...
					phihq(i, q) * whq(q);

				fE(i + 7) = fE(i + 7) + loading(xq, yq, 2) * ...
					phihq(i, q) * whq(q);
			end

        	fE(i) = 2 * area * fE(i);
        	fE(i + 7) = 2 * area * fE(i + 7);
		end
	
    	% Builds the matrices.
		A(dofgg, dofgg) = A(dofgg, dofgg) + KE;
    	B(dofp, dofgg) = BE;

		% Builds the RHS.
		b(dofgg) = b(dofgg) + fE;
		
		% Multipliers.
		areaVector(dofp) = area * ones(3, 1);
	end

	% Homogeneus DBCs.
	if ~(isempty(approximate)) % P2 approximation.
		% Strips bubble.
		dofs = verticesNumber + elementsNumber + edgesNumber;
		bubble = setdiff(1:1:dofs, verticesNumber + (1:1:elementsNumber));
			
		bubble = [bubble, verticesNumber + elementsNumber ...
			+ edgesNumber + bubble];

		A = A(bubble, bubble);
		B = B(:, bubble);
		b = b(bubble);
		
		% BCs on P2 approximation.
		nodesI = setdiff(1:1:verticesNumber, boundary);
		edgesI = setdiff(1:1:edgesNumber, boundedges);
		NL = [nodesI, verticesNumber + edgesI];
		NL = [NL, verticesNumber + edgesNumber + NL];

		uh = zeros(2 * (verticesNumber + edgesNumber), 1);

	else
		% Standard BCs.
		nodesI = setdiff(1:1:verticesNumber+elementsNumber, boundary);
		edgesI = setdiff(1:1:edgesNumber, boundedges);
		NL = [nodesI, verticesNumber + elementsNumber + edgesI];
		NL = [NL, verticesNumber + edgesNumber + elementsNumber + NL];

		uh = zeros(2 * (verticesNumber + edgesNumber + elementsNumber), 1);
	end

	% Matrices without DBCs.
	Ah = A(NL, NL); clear A;
	Bh = B(:, NL); clear B;
	fh = b(NL); clear b;

	% Sparse matrices.
	Ah = sparse(Ah);
	Bh = sparse(Bh);
	fh = [fh; zeros(3 * elementsNumber, 1)];

	% Global matrix.
	Kh = [Ah, Bh'; Bh, zeros(3 * elementsNumber)];
	
	% Multiplier conditions on homogeneus DBCs.
	N = length(NL);
	Kh = [Kh, [zeros(N,1); areaVector]; [zeros(1, N), areaVector', 0]];
	fh = [fh; 0];

	% Linear system solution.
	Kh = sparse(Kh); % Probably redundant.
	solh = Kh \ fh;

	% Extracts the two components of the velocity and
	% pressure.
	uh(NL) = solh(1:N);
	uh = reshape(uh, length(uh) / 2, 2);
	ph = solh(N+1:(N+(3*elementsNumber)));

	% Solution data.
	solution = struct;

	solution.uh = uh;
	solution.ph = ph;

	% Error evaluation.
	error = struct;

	l2Error = 0;
	h1Error = 0;
	l2ErrorP = 0;

	for index = 1:elementsNumber
		% Vertices information.
		vertex = [vertices(index, 1), vertices(index, 2), ...
			vertices(index, 3)];
		
		X = [xv(vertex(1)), xv(vertex(2)), xv(vertex(3))];
		Y = [yv(vertex(1)), yv(vertex(2)), yv(vertex(3))];

		% Jacobian.
		jacob = [X(2) - X(1), X(3) - X(1); Y(2) - Y(1), Y(3) - Y(1)];

		% Inverse Jacobian.
		invJacob = inv(jacob);

		% Transpose inverse Jacobian.
		invJacobT = invJacob';

		% Area.
		area = 0.5 * det(jacob);

		% Edges information.
		edge = [edges(index, 1), edges(index, 2), edges(index, 3)];

		% 'Local to global' arrays of dofs.
		if isempty(approximate)
			dofg = [vertex(1), vertex(2), vertex(3), ...
				verticesNumber + index, ...
				verticesNumber + elementsNumber + edge(1), ...
				verticesNumber + elementsNumber + edge(2), ...
				verticesNumber + elementsNumber + edge(3)];
		else
			dofg = [vertex(1), vertex(2), vertex(3), ...
				verticesNumber + edge(1), ...
				verticesNumber + edge(2), ...
				verticesNumber + edge(3)];
		end

		dofp = [0, elementsNumber, 2 * elementsNumber] + index;

		% Solutions at dofs.
		uT1 = uh(dofg, 1);
		uT2 = uh(dofg, 2);
		pT = ph(dofp);

		% Partials.
		l2Err = 0;
		h1Err = 0;
		l2ErrP = 0;

		% Quadrature.
		for q = 1:Nq
			% uh, first and second component, values.
			uh1 = 0;
			uh2 = 0;

			% uh, first and second component, gradient's values.
			uhg1 = [0; 0];
			uhg2 = [0; 0];
			
			% Velocity.
			j = 1;

			for i = 1:7
				if i == 4 & ~(isempty(approximate))
					continue;
				end

				uh1 = uh1 + uT1(j) * phihq(i, q);
				uh2 = uh2 + uT2(j) * phihq(i, q);

				uhg1 = uhg1 + uT1(j) * (invJacobT * ...
					[gphihqx(i, q); gphihqy(i, q)]);
				uhg2 = uhg2 + uT2(j) * (invJacobT * ...
					[gphihqx(i, q); gphihqy(i, q)]);

				j = j + 1;
			end

			% ph values.
			phv = pT(1) * pphihq(1, q) + ...
				pT(2) * pphihq(2, q) + ...
				pT(3) * pphihq(3, q);
			
			% Gauss node on the element.
			node = jacob * [xhq(q); yhq(q)] + [X(1); Y(1)];
			xq = node(1);
			yq = node(2);

			% Partial error.
			l2Err = l2Err + (exact(xq, yq, 1) - uh1) ^ 2 * whq(q) + ...
				(exact(xq, yq, 2) - uh2) ^ 2 * whq(q);

			h1Err = h1Err + norm(exact(xq, yq, 4) - uhg1) ^ 2 * whq(q) + ...
				norm(exact(xq, yq, 5) - uhg2) ^ 2 * whq(q);

			l2ErrP = l2ErrP + (exact(xq, yq, 3) - phv) ^ 2 * whq(q);
		end

		l2Err = l2Err * 2 * area;
		h1Err = h1Err * 2 * area;
		l2ErrP = l2ErrP * 2 * area;

		% Total error (squared).
		l2Error = l2Error + l2Err;
		h1Error = h1Error + h1Err;
		l2ErrorP = l2ErrorP + l2ErrP;
	end
	
	% Error data.
	error.l2 = sqrt(l2Error);
	error.h1 = sqrt(h1Error);
	error.l2P = sqrt(l2ErrorP);

	% Geometry data.
	geometry = struct;
	
	geometry.mesh = mesh;
	geometry.elements = elementsNumber;
	
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

	geometry.diameter = diameter;
	geometry.dofs = length(fh);

	% Plot, if requested.
	if ~(isempty(show))
		tiledlayout(2, 2);

		% Velocity plot on vertices dofs.
		% Numerical vs analytical.
		nexttile;

		quiver(xv, yv, uh(1:verticesNumber, 1), ...
			uh(1:verticesNumber, 2));

		hold on;
		title("Numerical velocity");
		hold off;
		
		nexttile;

		uex = zeros(verticesNumber, 1);
		uey = zeros(verticesNumber, 1);
	
		for i = 1:verticesNumber
			uex(i) = exact(xv(i), yv(i), 1);
			uey(i) = exact(xv(i), yv(i), 2);
		end
	
		quiver(xv, yv, uex, uey);

		hold on;
		title("Analytical velocity");
		hold off;

		% Pressure plot on centroids.
		nexttile;

		for k = 1:elementsNumber
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

		for k = 1:elementsNumber
			hold on;

			indexes = vertices(k, 1:3)';
			
			xc = sum(xv(indexes)) / 3;
			yc = sum(yv(indexes)) / 3;

			pressure = exact(xc, yc, 3) * ...
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