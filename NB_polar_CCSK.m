% q=64=> CCSK 0111011001011101011001110000010000101110000111011100100001101011
clear
pth1 = ( '.\related_functions');
addpath(pth1);
rng(0);
max_gen = 1e4; %MontCarlo simulation count
q=64;  %GF(q)
N = 64; %code length
K = 42; % information sequence length
SNRs_db=-7.5;
if q==64
    CCSK_seq_str = '0111011001011101011001110000010000101110000111011100100001101011'; %GF64
elseif q==128
    CCSK_seq_str = '10001010101000001101100110101101000001100011011001101100000101000011011110011100101010000111101001001110111111111000110111010100';%GF128
elseif q==256
    CCSK_seq_str = '0111010101010010010011011010111100100110001000010100100100011010101000010101011101101110000001100110011100001100100101100011100010011110011100101011110100111001100010110000101100000100100111011110000110000111010100010111111101111111011110111100000001101110';
elseif q==512
    CCSK_seq_str = '10100000110010000011000111111100011010100110101011010000101111100000011111010101101011000100010100000101011000001011011011011000101110100011001101000000000101111011101100100000011011111011111110001110101001100010000010111011000111010010110010000111100100101011110100011010001110000111100011110011001010100010101111001100011101101100100000110001000010000011010010001011100111111111010101000111000111100111100110011100100101111011100100100100111100101110011000110111101100111111011110110000000111110110010101001000';
end

pth = sprintf('.\\reliab_seq_diff_GF_diff_N_diff_SNR\\GF%d_ccsk_reliab_seq_L2H\\N%d\\', q, N);
fl_nm = sprintf('mat_N%d_GF%d_SNR%.3f.txt', N, q, SNRs_db);
data = read_reliability_seq([pth fl_nm]);




p=log2(q);

eta0 = CCSK_seq_str' - 48; %convert string to array of double
etaq = zeros(q,q); % table of all circularly rotated  eta0
for i=0:q-1
    etaq(:,i+1)=circshift(eta0,-i);
end
etaqm = 1-2*etaq; % modulation (1->-1 and 0->+1)
n=log2(N);

SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigma = sqrt(N0);
I1=bi2de(fliplr(de2bi((0:N-1)', n)));

reliab_seq = data.reliab_seq_prob+1;% +1 because if channel indexing startes from 0
%%


Hadamard = zeros(q,q); %Hadamard lookup table to simplify VN if its inputs already in freq domain
for i = 1:q
    Hadamard(i, i)=1;
    Hadamard(:,i) = fwht(Hadamard(:,i))*q;
end
count_correct_symb = zeros(N,1);

FE = 0;
FER = 0;

msg = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f\n", SNRs_db, 0, 0, 0); fprintf(msg)
for  gen_seq_cnt=1:max_gen
    info_seq = randi([0 q-1],N,1); %generate N random symbol to encode (N instead of K!!)
    u = info_seq;
    x=encode(u); %x is the encoded sequence to be transmitted
    % uu = inv_encode(x); inverse encoding
    frozen_symbols = u(reliab_seq(1:N-K));
    etax = zeros(q, N);
    for i = 1 : N
        etax(:,i)=etaqm(:,x(i)+1); % encoded sequence CCSK modulation
    end
    nse = sigma*randn(size(etax));
    y = etax + nse;
    % L  = (LLR_CCSK(y, etaq, q, N, sigma(i0)^2));
    L  = LLR_CCSK_FFT(y, eta0, N, sigma);
    %%
    decw = polar_nb_dec(L, Hadamard, reliab_seq, frozen_symbols, 1);
    FE = FE + ~isequal(decw,u);
    FER = FE/gen_seq_cnt;
    idx_corr = find(decw==u);
    count_correct_symb(idx_corr)=count_correct_symb(idx_corr)+1;
    if mod(gen_seq_cnt, 100)==0
        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f\n", SNRs_db, FE, gen_seq_cnt, FER);
        fprintf(msg)
    end
    

end

