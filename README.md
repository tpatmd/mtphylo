Parameters used
Minimum Number Of Sequences For A Conserved Position: 21
Minimum Number Of Sequences For A Flanking Position: 21
Maximum Number Of Contiguous Nonconserved Positions: 8
Minimum Length Of A Block: 5
Allowed Gap Positions: With Half
Use Similarity Matrices: Yes

translatorx.pl -i ${input_cds} -o ${input_cds.baseName} -p F -c 5 -t F -g "-b2=${params.gf} -b4=5 -b5=h"