
#Prepare after strelka2
for i in 909*; do bcftools view -fPASS $i/results/variants/somatic.snvs.vcf.gz | grep -v "#" | awk 'BEGIN{print("#CHROM","POS","REF","ALT","NAAF","NN","TAAF","TN","QSS_NT","SomaticEVS")}{ref=$4;alt=$5;refid=$4."U";altid=$5."U";nref=0;nalt=0;qssnt=gensub(/^.+;QSS_NT=([0-9]+);.*$/,"\\1","g",$8);somaticevs=gensub(/^.+SomaticEVS=/,"","g",$8);split($9,keys,":");for (i=1;i<=length(keys);i++){if(keys[i]==refid){nref=i}else if(keys[i]==altid){nalt=i}};split($10,normalds,":");split($11,tumords,":");normalrefs=normalds[nref];normalalts=normalds[nalt];tumorrefs=tumords[nref];tumoralts=tumords[nalt];split(normalrefs,normalref,",");split(normalalts,normalalt,",");split(tumorrefs,tumorref,",");split(tumoralts,tumoralt,",");nn=(normalref[1]+normalalt[1]);tn=(tumorref[1]+tumoralt[1]);if(nn>0){naf=normalalt[1]/nn}else{naf="NA"};if(tn>0){taf=tumoralt[1]/tn}else{taf="NA"};print $1,$2,$4,$5,naf,nn,taf,tn,qssnt,somaticevs}'>$i/${i}_freqs.tsv;done

#Prepare after multisnv
awk 'function nid(nc){if(nc=="A"){return 1}else if(nc=="C"){return 2}else if (nc=="G"){return 3}else{return 4}};function af(bcounts,ref,alt){split(bcounts,thisfundata,",");nref=thisfundata[nid(ref)];nalt=thisfundata[nid(alt)];if(nref+nalt>0){return(nalt/(nref+nalt))}else{return("NA")}};BEGIN{mcol=10}{if(/^#/){if(/^##/){}else{split($0,oheader);print $1,$2,$4,$5,$6,"SAMPLE","NAAF","NN","TAAF","TN","GQ","MTN","MTT","NSAMPLES"}}else{if($7=="PASS"){split($9,format,":");split($0,fields);for(i=1;i<=length(format);++i){if(format[i]=="GQ"){gqi=i}else if(format[i]=="SS"){somatici=i} else if(format[i]=="BCOUNT"){countsi=i} else if(format[i]=="DP"){depthi=i}};split($10,normaldata,":");nn=normaldata[depthi];naaf=af(normaldata[countsi],$4,$5);mtn=normaldata[somatici];ndist=0;;for(i=mcol+1;i<=NF;++i){split(fields[i],thisdata,":");if(thisdata[somatici]!=4 && thisdata[somatici] != mtn){ndist=ndist+1}};for(i=mcol+1;i<=NF;++i){split(fields[i],thisdata,":");taaf=af(thisdata[countsi],$4,$5);print $1,$2,$4,$5,$6,oheader[i],naaf,nn,taaf,thisdata[depthi],thisdata[gqi],mtn,thisdata[somatici],ndist}}}}' 909_multisnv.vcf > 909_multisnv.processed.tsv

#Prepare after lofreq
for i in 909*; do bcftools view -fPASS $i/*somatic_final_minus-dbsnp.snvs.vcf.gz | grep -v "#" | awk 'BEGIN{print("#CHROM","POS","REF","ALT","TAAF","TN","UQ")}{uq=gensub(/^.+;UQ=([0-9]+)$/,"\\1","g",$8);taaf=gensub(/^.+AF=([^;]+);.*$/,"\\1","g",$8);tn=gensub(/^DP=([^;]+);.*$/,"\\1","g",$8);print $1,$2,$4,$5,taaf,tn,uq}'>$i/${i}_freqs.tsv;done

gatk -T RealignerTargetCreator -I ZT-909-1/ZT-909-1.bam  -I ZT-909-2/ZT-909-2.bam -I ZT-909-3/ZT-909-3.bam -I ZT-909-5/ZT-909-5.bam -R $HUMAN_GENOME -o 909.intervals -known /home/dmalload/my_storage/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -known /home/dmalload/my_storage/1000G_phase1.indels.b37.vcf.gz

## Summarize coverage and duplication
echo "Sample, Good pairs (M pairs), Duplicate proportion, Library size (M molecules), %Exome at 5X, %Exome at 10X, %Exome at 20X, %Exome at 40X";for i in */coverage_exome; do sample=$(dirname $i);gpairs=$(grep "properly paired" $sample/qc_$sample/${sample}_flagstat.txt | awk '{print sprintf("%0.1f",$1/1000000)}');dupprop=$(cat $sample/${sample}_markduplicates.metrics | grep -v "#" | grep $sample | awk '{print sprintf("%.2f",$9)}');libsize=$(cat $sample/${sample}_markduplicates.metrics | grep -v "#" | grep $sample | awk '{print sprintf("%.1f",$10/1000000)}');breadth=$(awk -v sample="$sample" 'BEGIN{OFS=", "}$1 == sample {print $10, $9, $8, $7}' $sample/coverage_exome.sample_summary);echo $sample, $gpairs, $dupprop, $libsize, $breadth;done