% -----------------------------------------------------------------------------
% Name of file  : read_reliability_seq.m
% Description   : Reads reliability sequence data for non-binary polar codes
%                 from a formatted text file. The file contains ordered
%                 reliability indices, entropies, and error probabilities
%                 obtained from Monte Carlo simulations.
%
% Author        : Abdallah Abdallah (<abdallah.abdallah@univ-ubs.fr>)
% Organisation  : Lab-STICC, UMR 6284, Université Bretagne Sud
% Date          : October 13, 2025
% Licence       : CeCILL-B, see LICENSE file:
%                 https://cecill.info/licences/Licence_CeCILL-B_V1-en.html
%
% Note          : This work has been funded by the French ANR project MIOT
%                 under grant Projet-ANR-24-CE93-0017
%                 Web site: https://project.inria.fr/miot/
% -----------------------------------------------------------------------------
% History:
%   Creation    : 13/10/2025 - Initial version
% -----------------------------------------------------------------------------
% FUNCTION DESCRIPTION:
%   Reads a reliability-sequence text file used in non-binary polar coding.
%   The file contains several numeric vectors representing:
%       1. Channel order by entropy (least → most reliable)
%       2. Channel order by error probability (least → most reliable)
%       3. Channel entropies
%       4. Channel error probabilities
%       5. Number of Monte Carlo observations
%
% INPUT:
%   filename : string
%              Full path to the text file containing the reliability data.
%
% FILE FORMAT EXAMPLE:
%   0 1 2 4 3 5 6 7
%   0 1 2 4 3 5 6 7
%   0.00025 0.00000 0.00000 ...
%   0.00046 0.00000 0.00000 ...
%   N_obs: 64000
%
% OUTPUT (structure):
%   data.reliab_seq_entropy : sequence based on channel entropy
%   data.reliab_seq_prob    : sequence based on channel error probability
%   data.entropy            : vector of entropies per channel
%   data.errprob            : vector of error probabilities per channel
%   data.N_obs              : number of Monte Carlo simulations
% -----------------------------------------------------------------------------

function data = read_reliability_seq(filename)

    % --- Validate file existence ---
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % --- Read all non-empty lines ---
    lines = {};
    while ~feof(fid)
        L = strtrim(fgetl(fid));  % remove leading/trailing whitespace
        if isempty(L) || all(L == 0)
            continue;              % skip blank lines
        end
        lines{end+1} = L; %#ok<AGROW>
    end
    fclose(fid);

    % --- Validate file structure ---
    if numel(lines) < 5
        error('File "%s" does not contain enough valid data lines.', filename);
    end

    % --- Parse numeric data ---
    seq1    = str2num(lines{1}); %#ok<ST2NM> % Reliability by entropy
    seq2    = str2num(lines{2}); %#ok<ST2NM> % Reliability by error prob
    entropy = str2num(lines{3}); %#ok<ST2NM> % Entropy values
    errprob = str2num(lines{4}); %#ok<ST2NM> % Error probabilities
    N_obs   = sscanf(lines{5}, 'N_obs: %d'); % Monte Carlo count

    % --- Store in structured output ---
    data.reliab_seq_entropy = seq1;
    data.reliab_seq_prob    = seq2;
    data.entropy            = entropy;
    data.errprob            = errprob;
    data.N_obs              = N_obs;
end
