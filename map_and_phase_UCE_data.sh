#!/bin/bash
#author Tobias Hofmann (tobiashofmann@gmx.net)
sample=$(basename "$1");
reference2=$(echo "$2" | sed 's/.fasta//');
reference=$(echo "$reference2" | sed 's/.fst//');
refbase=$(basename "$reference");

mkdir ./$sample-results;
echo "starting mapping of $sample to $refbase";
/state/partition3/CLC-AssemblyCell/clc-assembly-cell-4.3.0-linux_64/clc_mapper -o ./$sample-results/$sample-$refbase-assembly.cas -d $2 -q $1/*READ1.fastq $1/*READ2.fastq -p fb ss 300 1000 --cpus 12;
echo "cas to bam convertion";
/state/partition3/CLC-AssemblyCell/clc-assembly-cell-4.3.0-linux_64/clc_cas_to_sam -a ./$sample-results/$sample-$refbase-assembly.cas -o ./$sample-results/$sample-$refbase.bam -f 33 -u;
echo "sorting";
/state/partition5/tobias/bin/samtools/samtools sort ./$sample-results/$sample-$refbase.bam ./$sample-results/$sample-$refbase.sorted;
echo "creating index file";
/state/partition5/tobias/bin/samtools/samtools index ./$sample-results/$sample-$refbase.sorted.bam;
echo "remove old bam";
rm ./$sample-results/$sample-$refbase.bam;
echo "**********$sample successfully mapped**********";
echo "====================================================";
echo "Starting phasing of $refbase for $sample";
echo "phasing bam files";
#creates two phased bam files called allele.0.bam and allele.1.bam
/state/partition5/tobias/bin/samtools/samtools phase -A -F -Q 20 -b ./$sample-results/$sample-$refbase.allele ./$sample-results/$sample-$refbase.sorted.bam;
echo "sorting phased bam files";
/state/partition5/tobias/bin/samtools/samtools sort ./$sample-results/$sample-$refbase.allele.0.bam ./$sample-results/$sample-$refbase.allele.0.sorted;
/state/partition5/tobias/bin/samtools/samtools sort ./$sample-results/$sample-$refbase.allele.1.bam ./$sample-results/$sample-$refbase.allele.1.sorted;
echo "mpileup and creating fq files";
/state/partition5/tobias/bin/samtools/samtools mpileup -u -f $2 ./$sample-results/$sample-$refbase.allele.0.sorted.bam | /state/partition5/tobias/bin/samtools/bcftools/bcftools view -cg - | /state/partition5/tobias/bin/samtools/bcftools/vcfutils.pl vcf2fq > ./$sample-results/$sample-$refbase.allele.0.fq;
echo "...";
/state/partition5/tobias/bin/samtools/samtools mpileup -u -f $2 ./$sample-results/$sample-$refbase.allele.1.sorted.bam | /state/partition5/tobias/bin/samtools/bcftools/bcftools view -cg - | /state/partition5/tobias/bin/samtools/bcftools/vcfutils.pl vcf2fq > ./$sample-results/$sample-$refbase.allele.1.fq;
echo "...";
/state/partition5/tobias/bin/samtools/samtools mpileup -u -f $2 ./$sample-results/$sample-$refbase.sorted.bam | /state/partition5/tobias/bin/samtools/bcftools/bcftools view -cg - | /state/partition5/tobias/bin/samtools/bcftools/vcfutils.pl vcf2fq > ./$sample-results/$sample-$refbase.original.fq;

echo "transforming fq into fasta";
seqtk seq -a ./$sample-results/$sample-$refbase.original.fq > ./$sample-results/$sample-$refbase.original.fasta;
seqtk seq -a ./$sample-results/$sample-$refbase.allele.1.fq > ./$sample-results/$sample-$refbase.allele.1.fasta;
seqtk seq -a ./$sample-results/$sample-$refbase.allele.0.fq > ./$sample-results/$sample-$refbase.allele.0.fasta;

rm $2.fai;

echo "**********completed $refbase for $sample**********";[
