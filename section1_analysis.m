
% define symbolic variables
syms A V r1 r2 r3 r4 k1 k2 

% non linear ODEs
dAdt = (r1 * V) + (r2 * A * V) + A * (r3 - (r4 * A));
dVdt = -(k1 * V) / (k2 + V);

fprintf('Calculating system equilibrium points \n\n');

% solve for equilibrium
sol = solve([dAdt == 0, dVdt == 0], [A, V]);

A_states = sol.A;
V_states = sol.V;

num_equilibria = length(A_states);

for idx = 1:num_equilibria
    fprintf('Equilibrium Point %d:\n', idx);
    fprintf('  A* = %s\n', char(A_states(idx)));
    fprintf('  V* = %s\n\n', char(V_states(idx)));
end

% use the average values
avg_r3 = 0.011606;
avg_r4 = 0.035875;

A_formula = A_states(2); 

% numerical verification
final_A = subs(A_formula, [r3, r4], [avg_r3, avg_r4]);

fprintf('Numerical Verification:\n');
fprintf('Long-Term Antibody Protection (A*) = %0.4f IU/mL\n', double(final_A));
