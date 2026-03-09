%This function implements a semi-supervised version of Belief Propagation

function [error, iterations] = semi_supervised_bp(A, na, cab, eta, target_error, max_iter)
% Created with the help of Gemini 3 Pro Preview
%
% Based on: "Asymptotic analysis of the stochastic block model..."
%             Decelle et al.
%
% Inputs:
%   A            : Adjacency matrix (NxN), sparse, symmetric (undirected).
%   na           : Prior probability for each group (1xq vector). Sum should be 1.
%   cab          : Affinity matrix (qxq).
%   target_error : Desired relative error.
%   max_iter     : Maximum number of iterations (e.g., 20).

    %% 1. Initialization and Pre-processing
    
    % Number of nodes
    N = size(A, 1);
    
    % Number of communities
    q = 2;
    
    % We assume N * na is integer-valued
    % Highest index in each community
    N1 = N * na(1);
    N2 = N;
    
    % Highest index of revealed label in each community
    R1 = round(N1 * eta);
    R2 = N1 + round((N2 - N1) * eta);
    
    % Ensure na is a row vector
    if iscolumn(na), na = na'; end
    
    % Create directed edge list from Adjacency matrix
    % (BP on undirected graph uses directed messages i->j and j->i)
    [rows, cols] = find(A); 
    E_directed = length(rows); % Total directed edges
    
    % Map edges to indices for vectorization
    % We need to know which index corresponds to the reverse edge j->i for every i->j
    % Encode edges as unique scalar IDs
    edge_ids = uint64(rows) + uint64(cols) * uint64(N); 
    rev_edge_ids = uint64(cols) + uint64(rows) * uint64(N);
    
    % Find indices
    [~, rev_map] = ismember(rev_edge_ids, edge_ids);
    
    % Indices of i -> j where label of i has been revealed
    disc_rows = rows <= R1;
    disc_ids1 = ismember(edge_ids, uint64(rows(disc_rows)) + uint64(cols(disc_rows)) * uint64(N));
    disc_rows = (rows > N1) & (rows <= N1 + R2);
    disc_ids2 = ismember(edge_ids, uint64(rows(disc_rows)) + uint64(cols(disc_rows)) * uint64(N));
    
    % Initialize Messages psi^{i->j} (randomly normalized)
    % Size: [Number of directed edges] x [q]
    psi = rand(E_directed, q);
    % Forcing disclosed information
    psi(disc_ids1, :) = [ones(sum(disc_ids1), 1) zeros(sum(disc_ids1), 1)];
    psi(disc_ids2, :) = [zeros(sum(disc_ids2), 1) ones(sum(disc_ids2), 1)];
    psi = psi ./ sum(psi, 2); % Normalize rows
    
    % Initialize External Field h (vector of size 1xq)
    % This represents the mean-field influence of non-edges (Eq. 27)
    % Starting with zeros or a small random perturbation
    h = zeros(1, q);                                                       
    
    %% 2. Belief Propagation Loop
    disp("Starting Semi-Supervised BP Inference");
    error = Inf;
    t = 1;
    while error > target_error && t < max_iter
        % --- Step A: Compute Interaction Terms (sum_{tk} c_{tk,ti} * psi_{k->i}) ---
        % This corresponds to the term inside the product in Eq. 22
        % psi_msg_in(e, ti) = Sum over tk of (c(tk, ti) * psi(e, tk))
        % Operation: Matrix multiplication (Edges x q) * (q x q)
        psi_msg_in = psi * cab; 
        
        % --- Step B: Aggregation at Nodes ---
        % Compute product of incoming messages for every node
        % To avoid underflow, we sum logarithms: log(Prod) = Sum(log(terms))
        log_psi_msg_in = log(psi_msg_in + 1e-20); % Add epsilon for stability
        
        node_log_sum = zeros(N, q);
        for g = 1:q
            % Accumulate incoming log-messages for group g at each target node
            % cols contains the target node indices for the edges
            node_log_sum(:, g) = accumarray(cols, log_psi_msg_in(:, g), [N, 1]);
        end
        
        % --- Step C: Compute Node Marginals (Eq. 28) ---
        % psi^i(ti) proportional to na(ti) * exp(-h(ti)) * Prod_incoming
        % We work in log domain first
        term_field = log(na + 1e-20) - h; % 1xq
        
        log_marginals = node_log_sum + term_field; % Broadcast add
        
        % Shift for numerical stability before exp (log-sum-exp trick)
        max_log = max(log_marginals, [], 2);
        marginals_unnorm = exp(log_marginals - max_log);
        
        % Normalize marginals
        Z_i = sum(marginals_unnorm, 2);
        node_marginals = marginals_unnorm ./ Z_i;
        
        % --- Step D: Update Field h (Eq. 27 / Line 11 logic) ---
        % h_ti = (1/N) * Sum_k Sum_tk c_{tk, ti} * psi^k(tk)
        % This is effectively the average interaction from all nodes
        avg_marginal = mean(node_marginals, 1); % 1xq
        h_new = avg_marginal * cab;
        
        % --- Step E: Update Messages (Cavity Method) (Eq. 26) ---
        % psi^{i->j} is proportional to NodeMarginal(i) / Interaction(j->i)
        % Working in log domain:
        % log(psi_new) = log(node_marginal_i) - log(interaction_from_j)
        % Note: node_log_sum(rows, :) gives the sum of ALL incoming to source i.
        % We subtract the specific incoming message from j (which is rev_map index)
        
        log_psi_new = (node_log_sum(rows, :) + term_field) - log_psi_msg_in(rev_map, :);
        
        % Normalize new messages
        max_val = max(log_psi_new, [], 2);
        psi_new = exp(log_psi_new - max_val);
        % Forcing disclosed information
        psi_new(disc_ids1, :) = [ones(sum(disc_ids1), 1) zeros(sum(disc_ids1), 1)];
        psi_new(disc_ids2, :) = [zeros(sum(disc_ids2), 1) ones(sum(disc_ids2), 1)];
        psi_new = psi_new ./ sum(psi_new, 2);
        
        % Updating psi and h
        psi = psi_new;
        h = h_new;
        
        % Group Assignments: argmax of marginals
        [~, group_assignment] = max(node_marginals, [], 2);
        
        % Error evaluation
        error = min([sum(abs(group_assignment - [ones(N1, 1); 2 * ones(N - N1, 1)])) sum(abs(group_assignment - [2 * ones(N1, 1); ones(N - N1, 1)]))]) / N;
        
        t = t + 1;
    end
    
    %% 3. Post-Processing
    
    % Number of iterations
    iterations = t - 1;

end