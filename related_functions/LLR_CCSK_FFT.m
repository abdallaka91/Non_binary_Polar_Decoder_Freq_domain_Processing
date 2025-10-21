% -----------------------------------------------------------------------------
% File Name     : LLR_CCSK_FFT.m
% Description   : Computes Log-Likelihood Ratios (LLRs) for CCSK modulation
%                 using FFT-based circular correlation in the frequency domain.
% Author        : Abdallah Abdallah (abdallah.abdallah@univ-ubs.fr)
% Organization  : Lab-STICC, UMR 6285, Université Bretagne Sud
% Date          : October 13, 2025
% License       : CeCILL-B, see LICENSE file:
%                 https://cecill.info/licences/Licence_CeCILL-B_V1-en.html
%
% Note          : Funded by the French ANR project MIOT (ANR-24-CE93-0017)
%                 Website: https://project.inria.fr/miot/
% -----------------------------------------------------------------------------
% History:
%   - Created: 13/10/2025
% -----------------------------------------------------------------------------
% Functional Description:
%
%   Computes CCSK Log-Likelihood Ratios (LLRs) in the frequency domain.
%   The received signal is correlated with all circular shifts of the
%   spreading sequence using FFT operations.
%
% -----------------------------------------------------------------------------
% INPUTS:
%   y_chan : [q × N] matrix of received CCSK symbols
%   PN     : [q × 1] base CCSK spreading sequence (binary/bipolar)
%   sigma  : scalar, noise standard deviation
%
% OUTPUT:
%   LLR    : [q × N] matrix of LLR values for all CCSK shifts
%
% -----------------------------------------------------------------------------
% ALGORITHM:
%   1.  Mirror PN sequence (for circular correlation) and apply FFT
%   2. Multiply PN and y_chan in the frequency domain
%   3. Apply inverse FFT to obtain correlations
%   4. Normalize and scale by (2/σ²)
% -----------------------------------------------------------------------------

function LLR = LLR_CCSK_FFT(y_chan, PN, sigma)

    % Step 1: Mirror PN sequence (for circular correlation) and apply FFT
    PN = [PN(1); PN(end:-1:2)];
    PN = fft(PN);

    % Step 2: FFT-based circular correlation
    LLR = ifft(PN .* fft(y_chan));

    % Step 3: Normalize and scale
    LLR = LLR - min(LLR, [], 1);
    LLR = LLR * (2 / sigma^2);

end
