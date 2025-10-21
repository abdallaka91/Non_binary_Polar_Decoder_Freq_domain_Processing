% -----------------------------------------------------------------------------
% File Name     : polar_nb_dec.m
% Description   : Non-binary polar decoder over GF(2^p) using frequency-domain
%                 belief propagation (hybrid CN–VN updates).
% Author        : Abdallah Abdallah (abdallah.abdallah@univ-ubs.fr)
% Organization  : Lab-STICC, UMR 6285, Université Bretagne Sud
% Date          : October 13, 2025
% License       : CeCILL-B, see LICENSE file:
%                 https://cecill.info/licences/Licence_CeCILL-B_V1-en.html
%
% Note          : This work has been funded by the French ANR project MIOT
%                 under grant ANR-24-CE93-0017
%                 Website: https://project.inria.fr/miot/
% -----------------------------------------------------------------------------
% History:
%   - Created: 13/10/2025
% -----------------------------------------------------------------------------
% Functional Description:
%
%   POLAR_NB_DEC performs soft-decision decoding of non-binary polar codes
%   defined over GF(2^p). The decoder alternates between:
%     • Check Node (CN) updates → executed in the frequency domain using FWHT.
%     • Variable Node (VN) updates → performed in either probability or frequency
%       domain depending on the layer and the data representation.
%
%   Each CN–VN kernel connects two branches (A, B) and produces two outputs:
%
%     CN_output ◄────⊕────── A
%                     │
%                     ▼
%     VN_output ◄────=────── B
%
%   - CN operation (⊕): GF(2^p) addition between A and B.
%   - VN operation (=): element-wise likelihood refinement.
%
%       • If A,B in probability domain:
%             permute(A ⊕ symbol) × B
%       • If A,B in frequency domain:
%             convert B → probability
%             keep A in frequency domain
%             multiply (Hadamard(:,symbol) ⊙ A)
%             convert result → probability
%             elementwise multiply with B
%
%   The “symbol” corresponds to the frozen or hard-decision symbol decided at
%   the CN output for the current kernel.
%
% -----------------------------------------------------------------------------
% INPUTS:
%   prb1           : [q × N] matrix, channel LLRs or probabilities per symbol
%   Hadamard       : [q × q] precomputed normalized FWHT basis
%   reliab_seq     : [1 × N] reliability sequence (least → most reliable)
%   frozen_symbols : [(N−K) × 1] frozen symbols (encoder side)
%   is_LLR         : boolean, true if prb1 contains LLRs, false if probabilities
%
% OUTPUT:
%   decw           : [N × 1] decoded codeword (symbols in GF(2^p))
%
% -----------------------------------------------------------------------------
% INTERNAL VARIABLES:
%   N, q, n   : code length, GF size, and log2(N)
%   V         : decision matrix per decoding layer
%   L         : 3D tensor of message probabilities (layer × node × q)
%   stat_v    : logical masks for processed clusters or nodes
%   is_freq   : indicator of frequency-domain representation per node
%   Int       : index map of input channels per layer i and node j
% -----------------------------------------------------------------------------
% Decoding tree node structure:
% Each node (i,j) at layer i splits into two child nodes at layer (i+1):
% one indexed (2*j−1) and the other (2*j). Each pair of sibling nodes at layer i
% shares a common parent node located at layer (i−1) with index (j+1)/2.
%
% The diagram below illustrates the hierarchical relationship between layers:
%   
%                        (i−1,(j+1)/2)
%                            │
%            ┌───────────────┴────────────────┐
%            │                                │
%            ▼                                ▼
%          (i,j)                            (i,j+1)
%      ┌─────┴──────┐                   ┌─────┴──────┐  
%      ▼            ▼                   ▼            ▼ 
%  (i+1,2*j−1)  (i+1,2*j)       (i+1,2*(j+1)−1)  (i+1,2*(j+1)) 
% -----------------------------------------------------------------------------

function decw = polar_nb_dec(prb1, Hadamard, reliab_seq, frozen_symbols, is_LLR)

N_K = length(frozen_symbols);
N   = size(prb1, 2);
q   = size(prb1, 1);
n   = log2(N);

% Convert from LLR to probability if needed
if is_LLR
    prb1 = exp(-prb1);
    prb1 = prb1 ./ sum(prb1, 1);
end

% Initialize decoded word and decision matrix
decw = nan(N,1);
V = nan(n+1, N);
V(n+1, reliab_seq(1:N_K)) = frozen_symbols;

% Propagate frozen values down through the decoding tree to avoid
% performing unecssary VN/CN operations
for l = n:-1:1
    for t = 0:N/2-1
        j = 2^(n-l);
        a = 2*t - mod(t,j) + 1;
        b = j + 2*t - mod(t,j) + 1;
        A = V(l+1, a); 
        B = V(l+1, b);
        if ~isnan(A+B)
            V(l, a) = bitxor(A, B);
            V(l, b) = B;
        end
    end
end

% Initialize likelihood containers
prob_3d = zeros(n+1, N, q);
prob_3d(1,:,:) = prb1';   % first layer = channel likelihoods

node_state = cell(n+1,1);
node_input_is_freq = cell(n,1);
node_inpts_indices = cell(n,1);
for i = 1:n+1
    nb_nodes = 2^(i-1);
    node_state{i,1} = false(1, nb_nodes); %at each layer i, there is 2^(i-1) nodes
    li = N/nb_nodes; %node input size
    if i <= n
        node_input_is_freq{i} = false(1, nb_nodes);
    end
    for j = 0:nb_nodes-1
        node_inpts_indices{i}(j+1,:) = j*li + (1:li); %node input channels indices
    end
end

% -------------------------------------------------------------------------
% Main decoding loop
% -------------------------------------------------------------------------
i = 1; j = 1; % i: layer index, j: node index
while i > 0
    if node_state{i}(j)  % if current node has already processed
        if mod(j,2) == 0 % if node is VN → backpropagate symbols to parent
            i = i - 1;
            j = j / 2;
            ii1 = node_inpts_indices{i}(j,:);
            temp0 = V(i+1,ii1);
            lg = length(ii1)/2;
            for k = 1:lg
                V(i,ii1(k))      = bitxor(temp0(k), temp0(k+lg));  % XOR backprop
                V(i,ii1(k+lg))   = temp0(k+lg);                    % equality path
            end
            node_state{i}(j) = true;
        else % if node is CN → go to its parent to process VN
            i = i - 1;
            j = (j + 1)/2;
        end

    elseif node_state{i+1}(2*j-1) % if left child processed → perform VN update
        ii1 = node_inpts_indices{i}(j,:);
        ii2 = node_inpts_indices{i+1}(2*j-1,:);
        lg  = length(ii1)/2;
        temp_s = V(i+1, ii2);  % symbol (frozen or hard-decision)

        % ---------- Variable Node (VN) Processing ----------
        if ~node_input_is_freq{i}(j) % VN in probability domain
            for k = 1:lg
                temp_a = reshape(prob_3d(i, ii1(k), :), [], 1);
                temp_b = reshape(prob_3d(i, ii1(k)+lg, :), [], 1);
                idx = bitxor(temp_s(k), 0:q-1) + 1;
                temp_c = temp_a(idx) .* temp_b;
                temp_c = temp_c / sum(temp_c);
                prob_3d(i+1, ii1(k)+lg, :) = temp_c;
            end
        else % VN in frequency domain
            for k = 1:lg
                temp_a = reshape(prob_3d(i, ii1(k), :), [], 1);
                temp_b = reshape(prob_3d(i, ii1(k)+lg, :), [], 1);
                temp_c = Hadamard(:, temp_s(k)+1) .* temp_a;
                temp_a = fwht(temp_a);
                temp_b = fwht(temp_b);
                temp_c = fwht(temp_c);
                prob_3d(i, ii1(k), :)       = temp_a;
                prob_3d(i, ii1(k)+lg, :)    = temp_b;
                temp_c = temp_c .* temp_b;
                temp_c = temp_c / sum(temp_c);
                prob_3d(i+1, ii1(k)+lg, :)  = temp_c;
            end
        end
        % ---------------------------------------------------

        node_input_is_freq{i}(j) = false; % VN output is always in probability domain
        i = i + 1; j = 2*j;
        if i <= n
            node_input_is_freq{i}(j) = false;
        end
        if i == n+1
            node_state{i}(j) = true;
            if isnan(V(n+1, j))
                [~, nn] = max(prob_3d(i,j,:));
                V(n+1,j) = nn - 1;
            end
        end

    else  % if left child not processed → perform CN update
        lg = length(node_inpts_indices{i}(j,:))/2;
        ii1 = node_inpts_indices{i}(j,:);
        if ~isnan(V(i+1, ii1))
            node_state{i}(j) = true;
        else
            if ~node_input_is_freq{i}(j)
                idx_all = [ii1(1:lg), ii1(1:lg)+lg];
                block = squeeze(prob_3d(i, idx_all, :)).';
                block = fwht(block) * q;
                prob_3d(i, idx_all, :) = permute(block, [3 2 1]);
                node_input_is_freq{i}(j) = true;
            end
            prob_3d(i+1, ii1(1:lg), :) = prob_3d(i, ii1(1:lg), :) .* prob_3d(i, ii1(1:lg)+lg, :);
            if i < n
                node_input_is_freq{i+1}(2*j-1) = true;
            end
            i = i + 1; j = 2*j - 1;

            if i == n+1
                node_state{i}(j) = true;
                if isnan(V(n+1, j))
                    temp = fwht(reshape(prob_3d(i,j,:), [], 1));
                    [~, nn] = max(temp);
                    V(n+1,j) = nn - 1;
                end
                i = i - 1; j = (j + 1)/2;
            end
        end
    end
end
decw = V(n+1,:)';
end
