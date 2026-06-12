% -----------------------------------------------------------------------------
% File Name     : polar_nb_encode_systematic.m
% Description   : Systematic non-binary polar encoder over GF(2^p)
%                 Uses non-systemtic encoder.m from A. Abdallah
%                 https://github.com/abdallaka91/Non_binary_Polar_Decoder_Freq_domain_Processing/blob/main/related_functions/encode.m
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
%   Algorithmic steps from 10.1109/TCOMM.2016.2574996:
%   1) v1 from u with zeros in frozen positions (Ac)
%   2) v2 = v1 · F^⊗m 
%   3) v3 = v2 with frozen positions (Ac) set to zero
%   4) x  = v3 · F^⊗m
%
%   Where F = [1 0; 1 1]
%
% -----------------------------------------------------------------------------
% INPUT:
%   u           : [K × 1] vector of GF(2^p) symbols (information sequence)
%   reliab_seq  : [1 × N] reliability sequence (least → most reliable)
%   N           : code length (power of 2)
%   K           : number of information symbols
%
% OUTPUT:
%   x           : [N × 1] vector of GF(2^p) symbols (systematic codeword)
%
% -----------------------------------------------------------------------------

function x = polar_nb_encode_systematic(u, reliab_seq, N, K)
    A = reliab_seq(N - K + 1:N);
    Ac = reliab_seq(1:N - K);

    u = reshape(u, [], 1);
    v1 = zeros(N, 1);
    v1(A) = u;

    v2 = encode(v1);

    v3 = v2;
    v3(Ac) = 0;

    x = encode(v3);
end