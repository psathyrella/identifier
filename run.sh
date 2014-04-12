#!/bin/bash
PATH=$PATH:../binaries
#hmmbuild v-genes.hmm ~/work/recombinator/data/ighv.fasta
hmmdir=~/work/recombinator/data/msa
for sto_file in `ls $hmmdir/*`; do
#    hmmbuild --informat stockholm --dna ighv5-a-01.hmm <(sed -e '1 s/^/# STOCKHOLM 1.0\n/' -e '$ s@$@\n//@' $sto_file)
    tmpfile=/tmp/`basename $sto_file`
    outfile=`echo $sto_file | sed 's/\.sto/\.hmm/'`
    sed -e '1 s/^/# STOCKHOLM 1.0\n/' -e '$ s@$@\n//@' $sto_file >$tmpfile
    hmmbuild --informat stockholm --dna $outfile $tmpfile
    rm -f $tmpfile
done
#cat $hmmdir/*.hmm >all-vdjs.hmm
#hmmpress all-vdjs.hmm
#hmmscan all-vdjs.hmm test-gene.fa
#7184821250938227179,5987052119601666025,IGHV3-23*01,IGHD4-17*01,IGHJ4*02_F,36,AT,TT,3,1,0,8,313AC-234CT-243CA-186TC-273GA-303TC-244AG-242GC-170CT-78TA-88GA-195GA-94AG-182AG-78AG-230CG-200GC-195AC-246AG-212CA-56AG-135GA-68AG-31TC,GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTCGGTACAGCCTGGGGGGTCCCTGAGGCTCTCCTGTGCGGCCTCTGGAGTCACCTTTAACAGCTGTGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGAAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGTACATACTACGCGGACCCCGTGAAGCGCCGCTTCACCATCTCAAGAGACAATTCCAAGAAGACGTTGTATCTCAGAGTGAACAGCCTGAGAGCCGAGGACACGACCGTATATTACTGTGCGAAATGACTACGGCGACTACTTGCCTACTGGGGCCAGGGA

#----------------------------------------------------------------------------------------
exit 1
for line in `cut -f3,4,5,14 -d, ~/work/recombinator/out.csv | sed '1 d'`; do
    v_gene=`echo $line | cut -f1 -d, | sed 's/\*/-/'`
    d_gene=`echo $line | cut -f2 -d, | sed 's/\*/-/'`
    j_gene=`echo $line | cut -f3 -d, | sed 's/\*/-/'`
    seq=`echo $line | cut -f4 -d,`
    tmpfile=/tmp/tmp-seq.fa
    echo "> $v_gene%$d_gene%$j_gene" >$tmpfile
    echo "$seq" >>$tmpfile
    hmmoutfile=/tmp/hmmout.txt
    hmmscan all-vdjs.hmm $tmpfile >$hmmoutfile
    inferred_v=`grep -m2 IGHV $hmmoutfile | grep -v Query | awk '{print $9}'`
    inferred_d=`grep -m2 IGHD $hmmoutfile | grep -v Query | awk '{print $9}'`
    inferred_j=`grep -m2 IGHJ $hmmoutfile | grep -v Query | awk '{print $9}'`
    printf "%15s%15s%15s" $v_gene $d_gene $j_gene
    if [ "$v_gene" != "$inferred_v" ]; then echo -n "   v wrong"; else echo -n "   v ok   "; fi
    if [ "$d_gene" != "$inferred_d" ]; then echo -n "   d wrong"; else echo -n "   d ok   "; fi
    if [ "$j_gene" != "$inferred_j" ]; then echo -n "   j wrong"; else echo -n "   j ok   "; fi
    echo ""
    rm -f $tmpfile $hmmoutfile
done

#nhmmer ighv5-a-01.hmm test-v-gene.fa
