import subprocess
import pandas as pd
import re
import os
from Bio.Seq import Seq

# Global variables:
mutation_size = 6
size_of_seq = 52
mut_num = 0 #max 100.  0 -> to get the common above 5%
mer_lenght = 9
counter = 0
list_of_HLA = ['HLA-A0201', 'HLA-C0701', 'HLA-A0101', 'HLA-C0702', 'HLA-C0401', 'HLA-A0301', 'HLA-B0702', 'HLA-A2402', 'HLA-B0801', 'HLA-C0602']

def most_common_mut(mut_num, primary_site):
    f_path = f"/home/alu/netlandes/MHCpan/frequent-mutations/{primary_site}.tsv"
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(f_path, delimiter='\t')
    # Filter rows where the third column does not start with "Stop Gained"
    filtered_df = df[~df.iloc[:, 2].str.startswith("Stop Gained")]

    if mut_num == 0:
        # Extract the fourth column from the DataFrame
        filtered_df[['Value1', 'Value2']] = filtered_df.iloc[:, 3].str.split('/', expand=True)
        filtered_df['Value2'] = filtered_df['Value2'].str.rsplit(',', n=1, expand=True)[1].str.replace('%', '').astype(float)
        filtered_df2 = filtered_df[filtered_df['Value2'] >= 5]              #!!!! ALL -> 1
        selected_rows = filtered_df2.iloc[:,0:4].values
    else:
        if mut_num > 100:
            mut_num = 100
        selected_rows = filtered_df.iloc[:mut_num, 0].values

    # Create a list of rows
    original_list = selected_rows.tolist()

    new_list = []

    for item in original_list:
        original_string = item[0]
        gene_mutation = item[2].split()[-2:]  # Get the last two words
        gene_name, mutation = gene_mutation
        precent = item[3].split(',')[1]                                      #!!!! ALL -> precent = item[3].split(',')[2]
        new_list.append([original_string, gene_name, mutation, precent])

    return new_list

def delete_files(folder_path):
    # iterate over all files in the folder
    for filename in os.listdir(folder_path):
        # construct the full file path
        file_path = os.path.join(folder_path, filename)
        
        # check if the file path is a regular file (not a folder or a symlink)
        if os.path.isfile(file_path):
            # delete the file
            os.remove(file_path)

def checking_input(mutation_input):
    # We want deletion and insertion mutations only:
    if 'del' in mutation_input:
        pattern = r'(\w+):g\.(\d+)([a-z]+)([A-Z]+)'
    elif 'ins' in mutation_input:
        pattern = r'(\w+):g\.(\d+)_\d+([a-z]+)([A-Z]+)'
    elif '>' in mutation_input:
        pattern = r'(\w+):g\.(\d+)(\w)>(\w)'
    else:
        return None
    return re.search(pattern, mutation_input)

def run_MHCpan(seq_list, patient_HLA, strand):
    count = 0
    a = ">Temp:"
    for seq in seq_list:
        count += 1
        mut_aa_seq = Seq(seq).translate()
        with open("/home/alu/netlandes/MHCpan/input.fsa", 'w') as f:
            f.write(a + '\n')
            f.write(str(mut_aa_seq))
            output_file = f"/home/alu/netlandes/MHCpan/results_common/result{count}.xls"  
        subprocess.run(f'"/private/common/Software/netMHC-4.0/Linux_x86_64/bin/netMHC" -f "/home/alu/netlandes/MHCpan/input.fsa" -tdir "/tmp" -version "/private/common/Software/netMHC-4.0/Linux_x86_64/data/version " -hlalist "/private/common/Software/netMHC-4.0/data/allelelist" -syn "/private/common/Software/netMHC-4.0/Linux_x86_64/data/synlists/%s.synlist" -thrfmt "/private/common/Software/netMHC-4.0/threshold/%s.thr" -rdir "/private/common/Software/netMHC-4.0/Linux_x86_64" -xlsfile {output_file} -xls -a {patient_HLA}', shell=True)

def create_results_file(tuples_of_SB, tuples_of_WB, hla):
    dir_files = os.listdir("/home/alu/netlandes/MHCpan/results_common/")
    i = 0
    for file in dir_files:
        path = "/home/alu/netlandes/MHCpan/results_common/" + file
        results = pd.read_csv(path, sep='\t', skiprows=1)
        new_results = results.loc[results['N_binders'] == 1]
        if not new_results.empty:
            for _, row in new_results.iterrows():
                if(row['Rank'] <= 0.5):
                    tuples_of_SB.append((i , row['Peptide'], row['Rank'], hla))
                else:
                    tuples_of_WB.append((i , row['Peptide'], row['Rank'], hla))
        i += 1


if __name__ == "__main__":
    with open ("/home/alu/netlandes/MHCpan/final_common.tsv", 'w') as final_file:
        final_file.write(f"Cancer type:\tGene-Name:\tAA-change:\tMutation:\tPrevalence-percentage\tBest-target:\tRank:\tHLA:\n")
        file_names  = os.listdir("/home/alu/netlandes/MHCpan/frequent-mutations")
        for primary_site in [os.path.splitext(file_name)[0] for file_name in file_names]:
            list_of_mutation = most_common_mut(mut_num, primary_site)
            for mut in list_of_mutation:
                mutation_input = mut[0]
                match = checking_input(mutation_input)
                if match == None: # A mutation that is not a deletion insertion Or Substitution
                    continue
                mut_chr = match.group(1)
                mut_pos = match.group(2)
                if '>' in mutation_input:
                    mut_type = "Substitution"
                else:
                    mut_type = match.group(3)
                mut_seq = match.group(4)

                # We want the size of the insertion or deletion to be no more than 6bp:
                if (len(mut_seq) > mutation_size):
                    print(f"The mutation should be a maximum {mutation_size}bp\n")
                    continue

                # Run a program that accepts as input a chromosome and location and returns the 60 bases around it:
                p = subprocess.run(f'sh "/home/alu/netlandes/MHCpan/find_seq_of_mutation.sh" {mut_chr} {mut_pos} {int(size_of_seq / 2)}', capture_output=True,text=True, shell=True)
                seq = p.stdout.replace("\n", "") 

                # We will run a program that will return the reading frame of the sequence:
                subprocess.run(f'sh "find.sh" {mut_chr} {mut_pos}', shell=True)
                output_file = "/home/alu/netlandes/MHCpan/intersects_output2.bed"
                if os.stat(output_file).st_size == 0:
                    print("The mutation is not in CDS\n")
                    continue
                temp = seq[::-1]
                isoforms_files = pd.read_csv(output_file, sep="\t", header=None)
                list_of_codons = []
                frame_init_pos = int(mut_pos) - int((size_of_seq/2) - 1)
                for index, row in isoforms_files.iterrows():
                    strand = row[4]
                    reading_frame = row[5]
                    if strand == "+":
                        cds_start_pos = row[1]
                        cds_end_pos = row[2]
                        count = 0
                        count2 =1
                    else:
                        count = 1
                        count2 =0
                        seq = temp
                        frame_init_pos = int(mut_pos) + int((size_of_seq/2) - 1)
                        cds_start_pos = row[2]
                        cds_end_pos = row[1]
                    end = len(seq)
                    flag = True
                    if abs(int(mut_pos) - cds_start_pos) < int((size_of_seq/2) - 1): 
                        start = int((size_of_seq/2) - 1) - abs(int(mut_pos) - cds_start_pos) + reading_frame + count
                        frame_init_pos = cds_start_pos + reading_frame 
                        flag = False
                    if abs(cds_end_pos - int(mut_pos)) < int(size_of_seq/2):
                        end = abs(cds_end_pos - frame_init_pos) 
                    if flag:
                        if (reading_frame == 0):
                            seq_frame = abs(frame_init_pos - cds_start_pos + count) % 3 
                        elif(reading_frame == 1):
                            seq_frame = abs(frame_init_pos - 1 - cds_start_pos) % 3 
                        else:
                            seq_frame = abs(frame_init_pos - 2 - cds_start_pos - count) % 3
                        if seq_frame == 0:
                            start = seq_frame 
                        if seq_frame == 1:
                            start = seq_frame + 1
                        if seq_frame == 2:
                            start = seq_frame - 1
                    codons = [seq[i:i+3] for i in range(start, end, 3)]
                    # Remove duplicates:
                    if codons not in list_of_codons:
                        list_of_codons.append(codons)
                # print(codons)

                # Creats mutatation in the seq for each isoform:
                mut_isoforms = []
                pos = int(size_of_seq/2) - start - count2 + count
                if strand == "-":
                    if mut_type == "ins":
                        mut_seq = mut_seq[::-1]
                    pos = pos - len(mut_seq)
 
                for codon in list_of_codons:
                    # Convert the list of codons to a DNA sequence
                    dna_seq = ''.join(codon)
                    if mut_type == "del":
                        # Remove the sequence from the DNA sequence:
                        mut_dna_seq = dna_seq[:pos] + dna_seq[pos + len(mut_seq):]
                    elif mut_type == "ins":
                        pos = int(size_of_seq/2) - start 
                        # Insert the sequence to the DNA sequence:
                        mut_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos:]
                    else:
                        mut_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos + 1:]
                    mut_isoforms.append(mut_dna_seq)

                ADAR_mut_seq_list = []
                for isoform in mut_isoforms:
                    # The netMHC tool accepts a peptide of at least 9mer in length
                    if len(isoform) < mer_lenght * 3:
                        continue
                    # If the mutation is on the negative strand, we will look for editing on the negative strand:
                    if strand == "-":
                        isoform = str(Seq(isoform).complement())
                    sequence_list = list(isoform)
                    # Iterate over the sequence and replace "A" with "G" when necessary
                    for i in range(1, len(sequence_list)):
                        sequence_list = list(isoform)
                        if sequence_list[i] == "A" and sequence_list[i - 1] != "G":
                            new_string = sequence_list
                            new_string[i] = "G"
                            ADAR_mut_seq_list.append("".join(new_string))

                # Remove duplicates from the adar_mut list:
                unique_ADAR_mut_list = list(set(ADAR_mut_seq_list))

                if unique_ADAR_mut_list:
                    sb_tuple = []
                    wb_tuple = []
                    for patient_HLA in list_of_HLA:
                        counter += 1
                        run_MHCpan(unique_ADAR_mut_list, patient_HLA, strand)
                        create_results_file(sb_tuple, wb_tuple, patient_HLA) 
                
                    delete_files('/home/alu/netlandes/MHCpan/results_common')                    

                    uniq_sb = list(set(sb_tuple))
                    uniq_wb = list(set(wb_tuple))

                    if not uniq_sb and not uniq_wb:
                        continue

                    best_ASO = ()
                    if uniq_wb:
                        for wb in uniq_wb:
                            wb_seq = Seq(unique_ADAR_mut_list[wb[0]])
                            protein = wb_seq.translate()
                            if "*" in str(protein):
                                protein = str(protein)
                                protein = protein.replace("*", "X")
                                protein = Seq(protein)
                            pos = protein.find(wb[1])
                            ASO_target = wb_seq[(pos * 3):((pos + len(wb[1])) * 3)]
                            best_ASO = (unique_ADAR_mut_list[wb[0]], ASO_target, wb[2], mutation_input, wb[3])
                            final_file.write(f"{primary_site}\t{mut[1]}\t{mut[2]}\t{best_ASO[3]}\t{mut[3]}\t{str(best_ASO[1])}\t{best_ASO[2]}\t{best_ASO[4]}\n")
                            
                    if uniq_sb:
                        for sb in uniq_sb:
                            sb_seq = Seq(unique_ADAR_mut_list[sb[0]])
                            protein = sb_seq.translate()
                            if "*" in str(protein):
                                protein = str(protein)
                                protein = protein.replace("*", "X")
                                protein = Seq(protein)
                            pos = protein.find(sb[1])
                            ASO_target = sb_seq[(pos * 3):((pos + len(sb[1])) * 3)]
                            best_ASO = (unique_ADAR_mut_list[sb[0]], ASO_target, sb[2], mutation_input, sb[3])
                            final_file.write(f"{primary_site}\t{mut[1]}\t{mut[2]}\t{best_ASO[3]}\t{mut[3]}\t{str(best_ASO[1])}\t{best_ASO[2]}\t{best_ASO[4]}\n")


