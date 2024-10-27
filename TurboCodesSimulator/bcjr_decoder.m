function LLR = bcjr_decoder(received_signal, trellis, prior_LLR)
    % BCJR Algorithm (MAP Decoder) for RSC code
    % received_signal: Noisy received signal
    % trellis: Trellis structure for the convolutional encoder
    % prior_LLR: A priori information (from extrinsic information of other decoder)
    
    N = length(received_signal);      % Length of signal
    num_states = length(trellis.nextStates); % Number of states in the trellis
    LLR = zeros(1, N);                % Pre-allocate LLRs (output)
    
    % Alpha and Beta (Forward and Backward state probabilities)
    alpha = zeros(num_states, N);     % Forward state metrics
    beta = zeros(num_states, N);      % Backward state metrics
    
    % Initialization
    alpha(:,1) = -Inf; alpha(1,1) = 0;  % Log(0) = -Inf, Log(1) = 0 for forward
    beta(:,N) = 0;                      % All backward state metrics initialized to 0
    
    % Forward Recursion (Alpha)
    for k = 2:N
        for state = 1:num_states
            % Calculate state transitions based on trellis
            for input = 0:1
                next_state = trellis.nextStates(state, input+1) + 1;  % Next state based on input
                branch_metric = received_signal(k-1) * trellis.outputs(state, input+1) + prior_LLR(k-1);  % Branch metric
                alpha(next_state, k) = logsumexp(alpha(next_state, k), alpha(state, k-1) + branch_metric);
            end
        end
    end
    
    % Backward Recursion (Beta)
    for k = N-1:-1:1
        for state = 1:num_states
            for input = 0:1
                next_state = trellis.nextStates(state, input+1) + 1;
                branch_metric = received_signal(k) * trellis.outputs(state, input+1) + prior_LLR(k);
                beta(state, k) = logsumexp(beta(state, k), beta(next_state, k+1) + branch_metric);
            end
        end
    end
    
    % Calculate LLRs based on Alpha and Beta
    for k = 1:N
        llr_numerator = -Inf;
        llr_denominator = -Inf;
        for state = 1:num_states
            for input = 0:1
                next_state = trellis.nextStates(state, input+1) + 1;
                branch_metric = received_signal(k) * trellis.outputs(state, input+1);
                alpha_beta_metric = alpha(state, k) + branch_metric + beta(next_state, k);
                if input == 1
                    llr_numerator = logsumexp(llr_numerator, alpha_beta_metric);  % Numerator for 1
                else
                    llr_denominator = logsumexp(llr_denominator, alpha_beta_metric);  % Denominator for 0
                end
            end
        end
        LLR(k) = llr_numerator - llr_denominator;  % LLR calculation (log domain)
    end
end

