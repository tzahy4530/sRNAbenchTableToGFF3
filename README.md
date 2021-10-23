# sRNAbenchTableToGFF3
this scripts using the prediction output of sRNAbench (novel.txt, novel451.txt) to create GFF3 file.

* Dependency:
  * Python 3.x
  * pandas

* Manual:

  `-i <path> : sRNAbench prediction output, like novel.txt/novel451.txt.`

  `-a <path> : additional input file.`
  
  `-o <path> : output path.`

* Example Run:

  `python sRNAbenchTableToGFF3 -i novel.txt -a novel451.txt -o sRNAbench_results.gff3`
