% Objective function for inverse problem
function J = obj_func_1Dheat_KL(c, data,alpha_mean, phi, lambda, f)
    alpha = alpha_mean + phi * (sqrt(lambda) .* c);
    model_sim = f(alpha);
    J = sum( (data - model_sim(end,:)).^2);

end