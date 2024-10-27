function decoded_bits = turbo_decoder(received_signal, trellis, iterations)
    % Turbo Decoder using BCJR Algorithm
    % received_signal: Noisy received signal (channel output)
    % trellis: Trellis structure for RSC encoder
    % iterations: Number of turbo iterations
    
    N = length(received_signal);  % Length of the received signal
    LLR = zeros(1, N);            % Log-likelihood ratios (soft information)
    decoded_bits = zeros(1, N);   % Pre-allocate decoded bits
    extrinsic_info = zeros(1, N); % Extrinsic information for iteration
    
    interleaver_map = randperm(N);
    deinterleaver_map = zeros(1, N);
    deinterleaver_map(interleaver_map) = 1:N;
    
    for iter = 1:iterations
        % First Decoder
        LLR_1 = bcjr_decoder(received_signal, trellis, extrinsic_info);  % Soft-input decoder
        extrinsic_info_1 = LLR_1 - extrinsic_info;  % Extrinsic information
        
        interleaved_extrinsic_info = extrinsic_info_1(interleaver_map);
        
        LLR_2 = bcjr_decoder(received_signal, trellis, interleaved_extrinsic_info);
        extrinsic_info_2 = LLR_2 - interleaved_extrinsic_info;  % Extrinsic information
        
        extrinsic_info = extrinsic_info_2(deinterleaver_map);
        
        decoded_bits = LLR_2 > 0;  
    end
end

