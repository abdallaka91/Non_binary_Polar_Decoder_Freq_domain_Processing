% -----------------------------------------------------------------------------
% File Name     : polar_nb_decode_systematic.m
% Description   : Systematic non-binary polar decoder over GF(2^p)
%                 Uses polar decoder from A. Abdallah
%                 https://github.com/abdallaka91/Non_binary_Polar_Decoder_Freq_domain_Processing/blob/main/related_functions/polar_nb_dec.m
%
% Author        : Sander Mikelsaar (mikelsaar.sander@univ-ubs.fr)
% Organization  : Lab-STICC, UMR 6285, Université Bretagne Sud
% Date          : June 11, 2026
% License       : CeCILL-B, see LICENSE file:
%                 https://cecill.info/licences/Licence_CeCILL-B_V1-en.html
%
% Note          : This work has been funded by the French ANR project MIOT
%                 under grant ANR-24-CE93-0017
%                 Website: https://project.inria.fr/miot/
%
% -----------------------------------------------------------------------------
% Functional Description:
%
%   1) decw = decode(y) 
%      (uses polar decoder to obtain decw)
%      (decw corresponds to v3 in the systematic encoder steps)
%   2) xHat = decw · F^⊗m 
%      (note that xHat != v1)
%   3) uHat = xHat(A) 
%      (information sequence from non-frozen positions A)
%
% -----------------------------------------------------------------------------
% INPUTS:
%   prb1           : [q × N] matrix, channel LLRs or probabilities per symbol
%   Hadamard       : [q × q] precomputed normalized FWHT basis
%   reliab_seq     : [1 × N] reliability sequence (least → most reliable)
%   N_K            : Number of frozen symbol
%   is_LLR         : boolean, true if prb1 contains LLRs, false if probabilities
%
% OUTPUT:
%   uHat           : [K × 1] decoded information sequence (symbols in GF(2^p))
%
% -----------------------------------------------------------------------------

function uHat = polar_nb_decode_systematic(prb1, Hadamard, reliab_seq, N_K, is_LLR)
    decw = polar_nb_dec(prb1, Hadamard, reliab_seq, N_K, is_LLR);
    decw = reshape(decw, [], 1);
    xHat = encode(decw);

    A  = reliab_seq(N_K + 1:end); % (Active positions)
    uHat = xHat(A);
end