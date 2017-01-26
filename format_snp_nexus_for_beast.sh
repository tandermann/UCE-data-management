cd topaza-uce-allele-snps
# the nexus file had to be edited a little to be able to use it in BEAST
#first I got the length of the alignment with this command:
length=$(fp.py --length snp.fasta | head -n 1)
#now create a list from 1 to x:
list=$(seq $length 1)
sed -i -e s/format.*\;/format\ datatype=standard\ symbols=\"012\"\ missing=?\;/g snp.nexus

sed -i -e '/missing=?;/a \
%		;\
' snp.nexus

for el in $list;
do sed -i -e "/missing=?;/a \\
%		$el	snp $el absent\/heterozygous\/homozygous\\
" snp.nexus
done

sed -i -e '/missing=?;/a \
%	CHARSTATELABELS\
' snp.nexus

sed -i '' 's/%//g' snp.nexus


rm snp.nexus-e
cd ..
