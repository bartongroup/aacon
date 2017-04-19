java -jar  compbio-conservation-1.1.jar -i=test/data/1000x3000DNA.aln.fa  -m=KABAT -o=kabat.txt -f=RESULT_NO_ALIGNMENT -d=kabat.d
java -jar  compbio-conservation-1.1.jar -i=test/data/large.aln.fa -m=JORES -o=jores.txt -f=RESULT_NO_ALIGNMENT -d=jores.d
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align  -m=KABAT,SCHNEIDER,Shenkin,gerstein,smerfs,valdar -o=schneider.txt -f=RESULT_NO_ALIGNMENT -d=schneider.d -n
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=SHENKIN -o=shenkin.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=GERSTEIN -o=gerstein.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=TAYLOR_GAPS -o=taylor.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=TAYLOR_NO_GAPS -o=taylor_nogaps.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=ZVELIBIL -o=zvelibil.txt -f=RESULT_NO_ALIGNMENT

java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=KARLIN -o=kalinin.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=ARMON -o=armon.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=THOMPSON -o=thompson.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=NOT_LANCET -o=not_lancet.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=MIRNY -o=mirny.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=WILLIAMSON -o=williamson.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=LANDGRAF -o=landgraf.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=SANDER -o=sander.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=VALDAR -o=valdar.txt -f=RESULT_NO_ALIGNMENT
java -jar  compbio-conservation-1.1.jar -i=test/data/TO1296.fasta.align -m=SMERFS -o=smerfs.txt -f=RESULT_NO_ALIGNMENT -d=smerfs.d -s=5,MID_SCORE,0.1
      
      