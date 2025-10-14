function decw = polar_nb_dec(prb1, Hadamard, reliab_seq, frozen_symbols, is_LLR)

%POLAR_NB_DEC  Non-binary polar decoder using probabilistic or LLR inputs.
%
%   decw = POLAR_NB_DEC(prb1, Hadamard, reliab_seq, frozen_symbols, is_LLR)
%
%   Performs non-binary polar decoding over GF(2^p) using belief propagation
%   in the frequency domain. The decoding alternates between Check Node (CN)
%   and Variable Node (VN) operations. CN updates are performed efficiently
%   using Fast Walsh–Hadamard Transforms (FWHT), while VN updates use
%   precomputed Hadamard-domain tables to avoid repeated FWHT computations.
%
%   INPUT ARGUMENTS:
%     1) prb1 : [q x N] matrix
%         - Each column corresponds to one symbol position.
%         - Each column is either:
%             • a probability vector (sum = 1), or
%             • a Log-Likelihood Ratio (LLR) vector.
%         - q = 2^p (the GF alphabet size).
%         - If LLRs are used, set is_LLR = true.
%
%     2) Hadamard : [q x q] matrix
%         - Contains the normalized FWHT of all GF(2^p) symbols.
%         - Each column corresponds to one symbol in the Hadamard domain.
%         - Example generation:
%
%             Hadamard = zeros(q,q);
%             for i = 1:q
%                 Hadamard(i,i) = 1;
%                 Hadamard(:,i) = fwht(Hadamard(:,i)) * q;
%             end
%
%     3) reliab_seq : [1 x N] vector
%         - Reliability sequence from least to most reliable bit channels.
%
%     4) frozen_symbols : [1 x (N−K)] vector
%         - Frozen symbols used during encoding, assigned as:
%             u(reliab_seq(1:N-K)) = frozen_symbols
%           where K is the number of information symbols (code rate = K/N).
%
%     5) is_LLR : logical flag
%         - true  → 'prb1' contains LLRs.
%         - false → 'prb1' contains probabilities (must be normalized).
%
%   OUTPUT:
%     decw : [1 x N] decoded codeword over GF(2^p).
%
%   INTERNAL DETAILS:
%     • Check Node (CN) operations are executed in the frequency domain using FWHT
%       to efficiently perform convolutions in GF(2^p).
%
%     • Variable Node (VN) updates reuse precomputed Hadamard-domain tables
%       (provided via the input 'Hadamard') to accelerate likelihood propagation
%       without requiring an explicit FWHT at every iteration.
%
%   EXAMPLE:
%       decw = polar_nb_dec(prb1, Hadamard, reliab_seq, frozen_symbols, true);
%
%   SEE ALSO:
%       FWHT, POLAR_ENC, POLAR_DEC, POLAR_NB_ENC
%
% -------------------------------------------------------------------------
% Author: Abdallah Abdallah
% Version: 1.1
% Description: Non-binary polar decoder (GF(2^p)) using hybrid CN–VN frequency-domain updates
% -------------------------------------------------------------------------

N_K = length(frozen_symbols);
N = size(prb1, 2);
q = size(prb1, 1);
n = log2(N);

if is_LLR
    prb1 = exp(-prb1);
    prb1 = prb1./sum(prb1,1);
end
decw = nan(N,1);
V = nan(n+1, N); %2D matrix contain the decoded symbol at each of the n+1 layer (layer 1 corrspond the the encoder input side)
V(n+1,reliab_seq(1:N_K)) = frozen_symbols;

for l=n:-1:1
    for t=0:N/2-1
        j = 2^(n-l);
        a=2*t-mod(t,j)+1;
        b=j+2*t-mod(t,j)+1;
        A = V(l+1, a); B = V(l+1, b);
        if ~isnan(A+B)
            V(l, a) = bitxor(A, B);
            V(l, b) = B;
        end
    end
end

L = zeros(n+1, N,q); %3D matrix of n+1 layers, the first layer contains the channel observation probabilities of each received symbol
L(1,:,:) = prb1'; %initialize 1st layer with channel observation input
stat_v = cell(n+1,1); %each layer contain 2^(l-1) cluster starting from layer l=0 corresponding to the decoder input side (channel side).
% stat_v{l}(cluster_idx is boolean 1 if it is treated, false if else.
for i = 1 : n+1
    i0 = i-1;
    stat_v{i,1} = false(1,  2^(i0)); % at the begining, all clusters aren't treated
end


is_freq = cell(n, 1); %is freq_is an indicator that tells if the data presented at the cluster inputs is in freq domain
Int = cell(n,1);%list on channel indexes of each cluster at each layer
for i = 1 : n+1

    mi = 2^(i-1);
    li = N/mi;
    if i<=n
        is_freq{i} = false(1, mi);
    end
    for j = 0:mi-1
        Int{i}(j+1,:) = j*li+1:j*li+1+li-1;
    end
end

i = 1;
j = 1;

while i>0
    if stat_v{i}(j)
        if mod(j,2)==0
            i=i-1;
            j=j/2;
            ii1 = Int{i}(j,:);
            temp0 = V(i+1,ii1);
            lg = length(ii1)/2;
            for i4=1:lg

                V(i,ii1(i4)) = bitxor(temp0(i4), temp0(i4+lg));
                V(i,ii1(i4+lg)) = temp0(i4+lg);
            end
            stat_v{i}(j) = true;
        else
            i = i-1;
            j=(j+1)/2;
        end
    elseif stat_v{i+1}(2*j-1)
        ii1=Int{i}(j,:);
        ii2=Int{i+1}(2*j-1,:);
        lg = length(Int{i}(j,:))/2;
        temp_s = V(i+1,ii2);
        temp_c = -ones(q,1);
        if ~is_freq{i}(j)
            for i4=1:lg
                temp_a =reshape(L(i,ii1(i4),:), [], 1);
                temp_b =reshape(L(i,ii1(i4)+lg,:), [], 1);

                for i5=0:q-1
                    temp_x = bitxor(temp_s(i4), i5);
                    temp_c(i5+1) = temp_a(temp_x)*temp_b(i5+1);
                end
                temp_c = temp_c/sum(temp_c);
                L(i+1,ii1(i4)+lg,:) = temp_c;
            end
        else
            for i4=1:lg

                temp_a =reshape(L(i,ii1(i4),:), [], 1);
                temp_b =reshape(L(i,ii1(i4)+lg,:), [], 1);
                temp_c = Hadamard(:, temp_s(i4)+1) .*temp_a;
                temp_a = fwht(temp_a);
                temp_b = fwht(temp_b);
                temp_c = fwht(temp_c);
                L(i,ii1(i4),:) = temp_a;
                L(i,ii1(i4)+lg,:) = temp_b;
                temp_c =  temp_c.*temp_b;
                temp_c = temp_c/sum(temp_c);
                L(i+1,ii1(i4)+lg,:) = temp_c;
            end
        end
        is_freq{i}(j) = false;
        i = i+1;
        j=2*j;
        if(i<=n)
            is_freq{i}(j) = false;
        end
        if i==n+1
            stat_v{i}(j)=true;
            if isnan(V(n+1, j))
                [mm,nn] = max(L(i,j,:));
                nn=nn-1;
                V(n+1,j) = nn;
            end
        end

    else
        lg = length(Int{i}(j,:))/2;
        ii1=Int{i}(j,:);
        if(~isnan(V(i+1, ii1)))
            stat_v{i}(j) = true;
        else
            if ~is_freq{i}(j)
                for i4=1:lg
                    L(i,ii1(i4),:) = fwht(reshape(L(i,ii1(i4),:), [], 1)) * q;
                    L(i,ii1(i4)+lg,:) =fwht(reshape(L(i,ii1(i4)+lg,:), [], 1)) * q;
                end
                is_freq{i}(j) = true;
            end
            for i4=1:lg

                temp_a = reshape(L(i,ii1(i4),:), [], 1);
                temp_b = reshape(L(i,ii1(i4)+lg,:), [], 1);

                temp_c = temp_a .*temp_b;
                L(i+1,ii1(i4),:) =  temp_c;
            end
            if(i<n)
                is_freq{i+1}(2*j-1) = true;
            end
            i=i+1;
            j=2*j-1;


            if i==n+1
                stat_v{i}(j)=true;
                if isnan(V(n+1, j))
                    temp = fwht(reshape(L(i,j,:), [], 1));
                    [mm,nn] = max(temp);
                    nn=nn-1;
                    V(n+1,j) = nn;
                end
                i=i-1;
                j=(j+1)/2;
            end
        end
    end
    decw = V(n+1,:)';
end
