1. Unpack all files from archive to the same folder. 

2. Copy sequence, fasta format, that you would like to map into this folder.

3. Convert sequence [file.fasta] from "fasta" format to "seq" format (one line) 
   by command:

   gawk -f fasta2seq.awk [file.fasta] > [file.seq]

4. Archive contains AA-TT dinucleotide matrix of human nucleosome, [AATT_human.mtr]. 
   You can use other mono or dinucleotide matrix with the same format.

5. Compile module mapping_CC.cpp by command:

   c++ mapping_CC.cpp -o mapping_CC

6. Change permission of script mapping_nucleosome_LinuxVer.s to executable
   by command:
   chmod 744 mapping_nucleosome_LinuxVer.s

7. Map nucleosome(s) on the sequence by command:

   mapping_nucleosome_LinuxVer.s [file_matrix.mtr] [file_sequence.seq] [file_output]

8. Your output [file_output] contains two columns: the first colume is coordinate and the 
   second colume is correlation coefficient (CC) of sequence to AATT dinucleotide matrix. 
   If sequence is long enough there may be more than one coordinate of nucleosome, 
   higher correlation is better.
