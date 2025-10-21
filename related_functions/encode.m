% -----------------------------------------------------------------------------
% File Name     : encode.m
% Description   : Recursive non-binary polar encoder over GF(2^p) 
%                 based on the Kronecker kernel F = [1 0; 1 1].
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
%   This function encodes a non-binary input vector u ∈ GF(2^p)^N
%   using the recursive Kronecker-based polar transform:
%
%       F = [1 0; 1 1]
%
%   The operation is performed symbol-wise using GF addition 
%   (implemented with bitwise XOR since GF(2^p) has characteristic 2):
%
%       x_a = u_a ⊕ u_b
%       x_b = u_b
%
%   The process is applied over n = log2(N) recursive stages.
%   The intermediate matrix `m1` stores partial results at each layer:
%     - Column 1     : encoder input (u)
%     - Column n + 1 : encoded output (x)
%
% -----------------------------------------------------------------------------
% INPUT:
%   u  : [N × 1] vector of GF(2^p) symbols (information + frozen symbols)
%
% OUTPUT:
%   x1 : [N × 1] encoded vector after full polar transformation
%
% -----------------------------------------------------------------------------

function x1 = encode(u)
N = length(u);
n = log2(N);

m1 = nan(N, n+1);    % Intermediate matrix storing layer states
m1(:,1) = u;         % Layer 1 = input symbols

% Recursive encoding: Kronecker kernel application
for l = 1:n
    for t = 0:N/2-1
        a = 2*t - mod(t, 2^(l-1)) + 1;
        b = 2^(l-1) + 2*t - mod(t, 2^(l-1)) + 1;
        tmp1 = [m1(a,l), m1(b,l)];
        m1(a,l+1) = bitxor(tmp1(1), tmp1(2));  % upper branch (⊕)
        m1(b,l+1) = tmp1(2);                   % lower branch (=)
    end
end

x1 = m1(:,end);      % Encoded vector (last layer output)
end
