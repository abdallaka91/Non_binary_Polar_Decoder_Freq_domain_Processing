% -----------------------------------------------------------------------------
% Name of file  : NB_polar_CCSK.m
% Description   : Monte Carlo simulation of Non-Binary Polar Codes (NB-PC)
%                 using CCSK modulation and frequency-domain decoding.
%                 This script simulates encoding, channel transmission, and
%                 decoding over GF(q) for a given SNR and code length.
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
% FUNCTIONAL DESCRIPTION:
%   This script performs Monte Carlo simulations of non-binary polar codes
%   using Cyclic Code Shift Keying (CCSK) modulation. The polar decoder is
%   implemented in the frequency domain for improved performance.
%
%   The simulation includes:
%       • Reliability sequence loading (from file)
%       • Polar encoding and CCSK modulation
%       • AWGN channel transmission
%       • Frequency-domain LLR computation
%       • Non-binary polar decoding (belief propagation)
%       • Frame Error Rate (FER) estimation
%
% -----------------------------------------------------------------------------
% INPUT PARAMETERS (set by user):
%   q         : Field order (GF(q)), must be a power of 2.
%   N         : Polar code length (number of encoded symbols).
%   K         : Information length (number of non-frozen symbols).
%   SNRs_db   : Channel Signal-to-Noise Ratio (in dB).
%   max_gen   : Number of Monte Carlo iterations.
%
% OUTPUT:
%   Console output with progressive FER estimation at given SNR.
% -----------------------------------------------------------------------------

%% --- INITIALIZATION --------------------------------------------------------
clear; clc;

if isunix
  addpath('./related_functions');  % Add path for related functions
else
  addpath('.\related_functions');  % Add path for related functions
end

rng(0);                          % Fix random seed for reproducibility

max_gen = 1e4;       % Monte Carlo simulation count
q = 64;              % GF(q)
N = 64;            % Code length
K = 42;             % Number of information symbols
N_K = N - K;         % Number of frozen symbols
SNRs_db = -7.5;       % SNR in dB

%% --- DEFINE CCSK MODULATION SEQUENCE ---------------------------------------
if q == 64
    CCSK_seq_str = '0111011001011101011001110000010000101110000111011100100001101011'; % GF(64)
elseif q == 128
    CCSK_seq_str = '10001010101000001101100110101101000001100011011001101100000101000011011110011100101010000111101001001110111111111000110111010100';
elseif q == 256
    CCSK_seq_str = '0111010101010010010011011010111100100110001000010100100100011010101000010101011101101110000001100110011100001100100101100011100010011110011100101011110100111001100010110000101100000100100111011110000110000111010100010111111101111111011110111100000001101110';
elseif q == 512
    CCSK_seq_str = '10100000110010000011000111111100011010100110101011010000101111100000011111010101101011000100010100000101011000001011011011011000101110100011001101000000000101111011101100100000011011111011111110001110101001100010000010111011000111010010110010000111100100101011110100011010001110000111100011110011001010100010101111001100011101101100100000110001000010000011010010001011100111111111010101000111000111100111100110011100100101111011100100100100111100101110011000110111101100111111011110110000000111110110010101001000';
else
    error('Unsupported value of q. Define CCSK sequence for this GF.');
end

%% --- LOAD RELIABILITY SEQUENCE ---------------------------------------------
if isunix
  pth = sprintf('.//reliab_seq_diff_GF_diff_N_diff_SNR//GF%d_ccsk_reliab_seq_L2H//N%d//', q, N);
else
  pth = sprintf('.//reliab_seq_diff_GF_diff_N_diff_SNR\\GF%d_ccsk_reliab_seq_L2H\\N%d\\', q, N);
end
fl_nm = sprintf('mat_N%d_GF%d_SNR%.3f.txt', N, q, SNRs_db);
data = read_reliability_seq([pth fl_nm]);

reliab_seq = data.reliab_seq_prob + 1;  % +1 if indexing starts from 0

%% --- BUILD MODULATION AND HADAMARD TABLES ----------------------------------
p = log2(q);
PN_bin = CCSK_seq_str' - 48;        % Convert binary string to double array
PN = 1 - 2 * PN_bin;                % Map 0→+1 and 1→-1

% Precompute circularly shifted CCSK symbols (modulation lookup table)
PN_mat = zeros(q, q);
for i = 0:q-1
    PN_mat(:, i+1) = circshift(PN, i);
end

% Precompute Hadamard transform table for VN frequency-domain updates
Hadamard = zeros(q, q);
for i = 1:q
    Hadamard(i, i) = 1;
    Hadamard(:, i) = fwht(Hadamard(:, i)) * q;
end

%% --- CHANNEL PARAMETERS ----------------------------------------------------
SNRs = 10.^(SNRs_db / 10);
N0 = 1 ./ SNRs;       % Noise power
sigma = sqrt(N0);     % Noise standard deviation

%% --- SIMULATION VARIABLES --------------------------------------------------
count_correct_symb = zeros(N, 1);
FE = 0;               % Frame Error counter
FER = 0;              % Frame Error Rate

msg = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f\n", SNRs_db, 0, 0, 0);
fprintf(msg);

%% --- MONTE CARLO SIMULATION LOOP -------------------------------------------
% ----- 1. initialize encoder input vector -----
u = zeros(N, 1);
for gen_seq_cnt = 1:max_gen
    % ----- 1. Information generation -----
    info_seq = randi([0 q-1], K, 1);
    u(reliab_seq(N_K+1:N)) = info_seq;

    % ----- 2. Encoding -----
    x = encode(u);

    % ----- 3. CCSK modulation -----
    x_mod = zeros(q, N);
    for i = 1:N
        x_mod(:, i) = PN_mat(:, x(i) + 1);
    end

    % ----- 4. AWGN channel -----
    nse = sigma * randn(size(x_mod));
    y = x_mod + nse;

    % ----- 5. LLR computation -----
    LLR = LLR_CCSK_FFT(y, PN_bin, sigma);
    % ----- 6. Decoding -----
    decw = polar_nb_dec(LLR, Hadamard, reliab_seq, N_K, true);

    % ----- 7. Error counting -----
    FE = FE + ~isequal(decw, u);
    FER = FE / gen_seq_cnt;

    % ----- 8. Progress display -----
    if mod(gen_seq_cnt, 10) == 0
        fprintf(repmat('\b', 1, length(char(msg))));
        msg = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f\n", SNRs_db, FE, gen_seq_cnt, FER);
        fprintf(msg);
    end
end

% -----------------------------------------------------------------------------
% END OF SCRIPT
% -----------------------------------------------------------------------------
