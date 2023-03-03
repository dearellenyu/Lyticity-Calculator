function lyticity_index = calculate_lyticity(sequence)

sequence_as_char = char(sequence);
% replace this sequence with your desired sequence
% sequence = 'TRSSRAGLQWPVGRVHRLLRK'

% segregate amino acids by hydrophilic vs hydrophobic
all_residues = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'];
hydrophobicity_value = [3.9,25.6,8.3,-0.9,-0.9,29.9,0.0,3.4,22.4,-1.1,24.2,16.3,0.5,9.7,0.5,3.9,0.5,3.9,14.4,32.9,29.5,15.4,-1.1];

% define hydrophilic residues
hydrophilic_residue = ['Q','N','S','R','D','H','K','Z','E'];

% make array of indices associating each residue of the peptide to its
% corresponding lyticity value
lyticity_array = [];
for i=1:strlength(sequence_as_char)
    amino_list_reference_number = find(ismember(all_residues,sequence_as_char(i)));
    lyticity_array(i) = amino_list_reference_number;
end

% make array of i+4 sums
i_plus_4_sums = [];
for i=1:(strlength(sequence_as_char)-5)
    if ismember(sequence_as_char(i),hydrophilic_residue) == 0
        if ismember(sequence_as_char(i+4),hydrophilic_residue) == 0
            i_plus_4_sums(i) = hydrophobicity_value(lyticity_array(i)) + hydrophobicity_value(lyticity_array(i+4));
        end
    end
end

% make array of i+3 sums
i_plus_3_sums = [];
for i=1:(strlength(sequence_as_char)-4)
    if ismember(sequence_as_char(i),hydrophilic_residue) == 0
        if ismember(sequence_as_char(i+3),hydrophilic_residue) == 0
            i_plus_3_sums(i) = hydrophobicity_value(lyticity_array(i))+hydrophobicity_value(lyticity_array(i+3));
        end
    end
end

% sum everything together to get lyticity value
lyticity_index = sum(i_plus_4_sums) + sum(i_plus_3_sums);

end
