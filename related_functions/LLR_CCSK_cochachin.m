% -----------------------------------------------------------------------------
% File Name     : LLR_CCSK_cochachin.m
% Description   : Computes the vector of Log-Likelihood Ratios (LLRs)
%                 for each CCSK-modulated GF(q) symbol using the
%                 Cochanchin analytical formulation under a BI-AWGN channel.
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
%   Implements the LLR definition proposed in:
%   Cochanchin V. et al., “Definition of the Log-Likelihood Ratio for
%   non-binary CCSK modulation over BI-AWGN,” IEEE Trans. on …
%   DOI: 10.1109/XXXX.9594140
%
%   For each received vector y(:,i), the function computes:
%
%       L_i(xi) = Σ_k [ (2 / σ²) * y_i(k) * (HD(k) − PN_xi(k)) ]
%
%   where:
%       - HD is the hard decision (0 if positive, 1 else) of the received symbol,
%       - PN_xi is the CCSK spreading sequence for symbol i,
%       - σ² is the noise variance.
%
%   The output L is normalized by subtracting its minimum element to
%   guarantee numerical stability (ensuring L ≥ 0).
%
% -----------------------------------------------------------------------------
% INPUTS:
%   y_chan : [q × N] matrix of received CCSK channel samples
%   PN     : [q × 1] fundamental binary PN (CCSK base) sequence
%   sigma  : noise standard deviation (√(N0))
%
% OUTPUT:
%   L      : [q × N] matrix of computed LLRs for each GF(q) symbol
%
% -----------------------------------------------------------------------------


function LLR = LLR_CCSK_cochachin(y_chan, PN, sigma)

% Precompute constants
N = size(y_chan, 2);           % Number of symbols (frames)
q = size(y_chan, 1);           % GF(q) order
Yscaled = (2 / sigma^2) * y_chan; % Scaled received samples
LLR = zeros(q, N);             % Output matrix
HD = y_chan < 0;               % Hard decision matrix η̂ (0 or 1)

% Main LLR computation per CCSK symbol
for i = 1:N
    for j = 1:q
        eta = circshift(PN(:), j-1) - HD(:,i);         % η_xi - η̂
        LLR(j, i) = Yscaled(:,i).' * eta;              % Eq. (3)
    end
    % Normalize to ensure non-negative LLRs
    LLR(:, i) = LLR(:, i) - min(LLR(:, i));
end
end
