# Author: Lei Wang
# Date: 10/12/2018
# EECS 730

1. After you are in this folder, you can type
> python main <filename1> <filename2>
to compile the program

2. I have included the human and mouse fasta files in the folder, you need to copy and paste files that you wish to test additionally within the folder. I also included the BLOSUM62.txt for reading matrix scores.

3. The result should display both the local and global alignments after applied affine gap penalty, it contains the alignment sequence and the corresponding scores.

4. I have commented out the intermediate matrixes such as the matrix after filled with scores and the matrix for trace backs. If you wish to check those information, you can simple uncomment line 106-110 and line 184 - 188. Then you can recompile the program and you will see the intermediate matrixes.

