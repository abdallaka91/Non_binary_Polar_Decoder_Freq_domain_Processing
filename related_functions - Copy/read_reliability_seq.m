function data = read_reliability_seq(filename)
    % READ_RELIABILITY_SEQ reads reliability sequences from a text file.
    %
    % File format:
    %  1. least-to-most reliable channels (entropies)
    %  2. least-to-most reliable channels (error probabilities)
    %  3. channels entropies
    %  4. channels error probabilities
    %  5. number of Monte Carlo simulations (N_obs)
    %
    % OUTPUT (structure):
    %  data.relab_seq_entropy - sequence based on channel entropy
    %  data.relab_seq_prob    - sequence based on channel error probability
    %  data.entropy           - entropies
    %  data.errprob           - error probabilities
    %  data.N_obs             - Monte Carlo simulations count

    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    lines = {};
    % Read all non-empty lines
    while ~feof(fid)
        L = strtrim(fgetl(fid));    % remove leading/trailing spaces
        if isempty(L) || all(L == 0)
            continue;               % skip blank lines
        end
        lines{end+1} = L; %#ok<AGROW>
    end
    fclose(fid);

    % Now parse based on expected structure
    if numel(lines) < 5
        error('File does not contain enough valid lines.');
    end

    seq1    = str2num(lines{1}); %#ok<ST2NM>
    seq2    = str2num(lines{2}); %#ok<ST2NM>
    entropy = str2num(lines{3}); %#ok<ST2NM>
    errprob = str2num(lines{4}); %#ok<ST2NM>
    N_obs   = sscanf(lines{5}, 'N_obs: %d');

    % Store in struct
    data.relab_seq_entropy = seq1;
    data.relab_seq_prob    = seq2;
    data.entropy           = entropy;
    data.errprob           = errprob;
    data.N_obs             = N_obs;
end
