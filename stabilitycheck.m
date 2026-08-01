% define symbolic variables
syms A V r1 r2 r3 r4 k1 k2

% non linear ODEs
f1 = (r1 * V) + (r2 * A * V) + A * (r3 - (r4 * A));
f2 = -(k1 * V) / (k2 + V);

fprintf('Calculating symbolic Jacobian matrix \n');

% compute the jacobian matrix for stability analysis
J = jacobian([f1; f2], [A, V]);
disp(J);

fprintf('Evaluating Jacobian at long-term immune state \n');

% substitute the long-term immune state (r3/r4, 0) into the jacobian
J_immune = subs(J, [A, V], [r3/r4, 0]);
disp(J_immune);

fprintf('Solving for symbolic eigenvalues... \n');

% find the eigenvalues to prove system stability
lambda = eig(J_immune);
disp(lambda);

% use the average parameter values
avg_r3 = 0.011606;
avg_r4 = 0.035875;
avg_k1 = 10.0;
avg_k2 = 50.0;

% numerical verification of the eigenvalues
numerical_lambda = subs(lambda, [r3, r4, k1, k2], [avg_r3, avg_r4, avg_k1, avg_k2]);

fprintf('Calculated Numerical Eigenvalues:\n');
fprintf('  Lambda 1 = %0.4f\n', double(numerical_lambda(1)));
fprintf('  Lambda 2 = %0.4f\n', double(numerical_lambda(2)));