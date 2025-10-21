% -----------------------------------------------------------------------------
% File Name     : inv_encode.m
% Description   : Inverse of the recursive non-binary polar encoder 
%                 over GF(2^p) using the kernel F = [1 0; 1 1].
%
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
% Functional Description:
%
%   Reconstructs the input vector u ∈ GF(2^p)^N from its encoded version x1
%   by recursively applying the **inverse polar transform**:
%
%       u_a = x_a ⊕ x_b
%       u_b = x_b
%
%   The process is identical to the encoding structure but executed 
%   in reverse order (from layer n down to 1).
%
%   The intermediate matrix `m1` stores partial decoding states:
%     - Column n + 1 : encoded symbols (decoder input)
%     - Column 1     : recovered original symbols
%
% -----------------------------------------------------------------------------
% INPUT:
%   x1 : [N × 1] vector of encoded GF(2^p) symbols
%
% OUTPUT:
%   u  : [N × 1] recovered vector of information + frozen symbols
%
% -----------------------------------------------------------------------------

function u = inv_encode(x1)
N = length(x1);
n = log2(N);

m1 = nan(N, n+1);
m1(:,n+1) = x1;  % Start from encoded vector

% Reverse recursion (inverse kernel application)
for l = n:-1:1
    for t = 0:N/2-1
        a = 2*t - mod(t, 2^(l-1)) + 1;
        b = 2^(l-1) + 2*t - mod(t, 2^(l-1)) + 1;
        tmp1 = [m1(a,l+1), m1(b,l+1)];
        m1(a,l) = bitxor(tmp1(1), tmp1(2));   % undo ⊕ operation
        m1(b,l) = tmp1(2);                    % propagate lower branch
    end
end

u = m1(:,1);  % Reconstructed input
end
