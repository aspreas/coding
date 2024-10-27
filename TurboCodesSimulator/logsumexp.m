function result = logsumexp(x1, x2)
    % Computes log(exp(x1) + exp(x2)) in a numerically stable way
    if x1 == -Inf
        result = x2;
    elseif x2 == -Inf
        result = x1;
    else
        max_val = max(x1, x2);
        result = max_val + log(exp(x1 - max_val) + exp(x2 - max_val));
    end
end
