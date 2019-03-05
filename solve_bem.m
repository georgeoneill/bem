function bem = solve_bem(surf,conductivity)

bem = struct;
bem.sigma = conductivity; %Condutvitity from outside to inside.
bem = add_gamma_multipliers(bem);
coeffs = fwd_bem_lin_pot_coeff(surf);
disp('Inverting Solution')
bem.solution = fwd_bem_multi_solution(coeffs, bem.gamma, [surf.np]);
if length(surf) == 3 &&  bem.sigma(2)/bem.sigma(3) < 0.1
    % Take the innermost compartment, recalculate coefficients and modify
    % the original solution accordingly.
    disp('IP Approach required')
    coeffs = fwd_bem_lin_pot_coeff(surf(3));
    disp('Inverting inner most shell, again')
    IP = fwd_bem_multi_solution(coeffs, [], surf(3).np);
    disp('Modifying solution')
    bem.solution = fwd_bem_ip_modify_solution(bem.solution,IP,bem.sigma(2)/bem.sigma(3),[surf.np]);
end
disp('BEM Solved')