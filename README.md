randomrnasnp
============

Add random snps and check for RNA structural changes

To use 

./pythonsnp24.py -f <fasta_file> -p <p_value cut off for RNAsnp1.1> -c <cutoff> -i <iterations> -q <yes/no>

where,

-f / --fasta : path to your fasta file containing your lncRNA sequences.
-p / --pvalue : p-value for RNAsnp 1.1.
-i / --iterations : iterations, number of times a single transcript is subjected to snp changes.
-c / --cutoff : cut off. number of transcript to have p-value lower than listed. ie, -p 0.5 -c 950 -i 1000 will suggest that 
                the program will stop when 950 iterations of the 1000 iterations have p-value lower than 0.5
-q / --quick : yes or no. 
                yes = Faster, stop the iterations when number of iterations with p-value above the threshold is
                higher than the (iterations - cutoff,)
                no = will continue the iteration. 
                




