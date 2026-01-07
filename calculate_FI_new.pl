# calculate Frequency of Impact for the RE (FI)
# for INV and INS, the impact is always 0.5, since the RE dosage is most probably not changed
# if in the same sample one RE was impacted in different levels, the more severe levels should be also counted for the less severe levels: e.g. if in one sample, RE1 has 100%impact for 1 allele, then it should also have 50% and <50% impact for 1 allele, unless it already have >=1 record of other impact levels

use warnings;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require SV_SUB;

$sv_file = "data/kg/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.hg19.vcf";
$re_anno = "data/reg_anno/combined_regulatory_elements_for_SVint.bed"; #sorted

undef %to_use; undef %sample2gt;
open F, $sv_file || die;
while(<F>){
	chomp; $line = $_;
	next if $line =~ /^##/;
	if($line =~ /^#/){
		$sample_no = scalar(split /\t/, $line) - 9; next;
	}
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  samples...
	($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = split /\t/, $_;
	next if $filter ne "PASS";
	$chr="chr$chr"; $stt = $pos - 1;
	next if $info !~ /END=/; next if $info !~ /SVTYPE=/;
	if($info =~ /END=(\d+);/){
		$end = $1;
	}else{
		$info =~ /SVLEN=(\d+);/; $end = $stt + $1;
	}
	$info =~ /SVTYPE=(.*)/; $tp=$1; # ALU;CNV;DEL;DEL_ALU;DEL_HERV;DEL_LINE1;DEL_SVA;DUP;INS;INV;LINE1;SVA
		$var = "$chr\t$stt\t$end";
# DEL
	if($tp =~ /^DEL/){
		$to_use{"DEL"}{$var}=1;
	}elsif($tp =~ /DUP/){
# DUP / CNV
		$to_use{"DUP"}{$var}=1;
	}elsif($tp =~ /CNV/){
#CNV
		$to_use{"CNV"}{$var}=1;
	}elsif($tp =~ /INV/){
# INV
		$to_use{"INV"}{$var}=1;
	}else{
# INS, ALU, LINE1, SVA
$var = "$chr\t$stt\t$pos";
		$to_use{"INS"}{"$chr\t$stt\t$pos"}=1; #$line;
	}
	@tmp = split /\t/, $line;
	for $i(9 .. $#tmp){
		$gt = $tmp[$i];
		$sum_gt = 0;
		if($gt ne "."){
			for $tt(split /\|/, $gt){
				$sum_gt = $sum_gt + 1 if $tt >= 1;
			}
		}
		next if $sum_gt == 0; # do not record 
		$sum_gt =2 if $sum_gt > 2;
		$sample = $i - 8;		
		$sample2gt{$var}{$sample} = $sum_gt;
	}
}
close F || die;

#open T, ">tmp3.txt";
#for $tp(sort keys %to_use){
 #       for $var(sort keys %{$to_use{$tp}}){
  #              print T "$tp\t$var\t$to_use{$tp}{$var}\n";
   #     }
#}
#close T;

#die;

&SV_SUB::write_bed2(\%to_use, "$sv_file\_temp");
$sv2re = SV_SUB::overlap("$sv_file\_temp.sorted.bed", $re_anno);
%sv2re = %{$sv2re};
unlink "$sv_file\_temp.sorted.bed";
unlink "$sv_file\_temp.bed";

# $..{$tp}{$sample_id}{$re}; for each sample, the same type of impact on the same RE will only be count once; for the same sample, different impact will all be taken into consideration
open O, ">$sv_file\_based_FI.txt" || die;
print O "SVTYPE\tCHR\tSTT\tEND\tFI(100% Affected)\tFI(50-100% Affected)\tFI(0-50% Affected)\n";
undef %all_1; # $all_1{$re_loc}{$sample}
undef %all_2; # $all_2{$re_loc}{$sample}
undef %all_3; # $all_3{$re_loc}{$sample}
for $tp(keys %to_use){
	print "$tp:\n";
	undef %count; # $count{$re_loc}{$im}
	for $i(1 .. $sample_no){
		undef %re2count;
		print "$i\n";
		for $var(keys %{$to_use{$tp}}){
#$line = $to_use{$tp}{$var};
#@tmp = split /\t/, $line;
#$col = 9 + $i -1; $gt = $tmp[$col];
			next if (!exists $sample2gt{$var}{$i});
			undef $sum_gt;			
			$sum_gt = $sample2gt{$var}{$i};
			for $re_loc(keys %{$sv2re{$var}}){
				$ratio = $sv2re{$var}{$re_loc};
				if($ratio >= 1){
					$im = "set3";
				}elsif($ratio < 0.5){
					$im = "set1";
				}else{
					$im = "set2";
				}
				if($tp eq "INS"){
					$im = "set2";
				}
# 100% impact (set3) will also be counted for 0-0.5 (set1) and 0.5-1 (set2); 0.5-1 (set2) impact will also be counted for 0-0.5 (set1);
				$re2count{$re_loc}{"set3"} = 0 if(!exists $re2count{$re_loc}{$im});
				$re2count{$re_loc}{"set2"} = 0 if(!exists $re2count{$re_loc}{$im});
				$re2count{$re_loc}{"set1"} = 0 if(!exists $re2count{$re_loc}{$im});
# theoretically each sample can only have <=2 disrupted alleles; will cut at 2 in the next block of scripts
				if($im eq "set3"){
					$re2count{$re_loc}{"set3"} += $sum_gt;
					$re2count{$re_loc}{"set2"} += $sum_gt;
					$re2count{$re_loc}{"set1"} += $sum_gt;
				}elsif($im eq "set2"){
					$re2count{$re_loc}{"set2"} += $sum_gt; $re2count{$re_loc}{"set1"} += $sum_gt;
				}else{
					$re2count{$re_loc}{"set1"} += $sum_gt;
				}
			}
		}
# for each case, the same re_loc with each impact should at most be counted twice (diplotype) 
		for $re_loc(keys %re2count){
			$count{$re_loc}{"set1"} = 0 if (!exists $count{$re_loc}{"set1"});
			$count{$re_loc}{"set2"} = 0 if (!exists $count{$re_loc}{"set2"});
			$count{$re_loc}{"set3"} = 0 if (!exists $count{$re_loc}{"set3"});
			$all_1{$re_loc}{$i} = 0 if (!exists $all_1{$re_loc}{$i});			
			$all_2{$re_loc}{$i} = 0 if (!exists $all_2{$re_loc}{$i});              
			$all_3{$re_loc}{$i} = 0 if (!exists $all_3{$re_loc}{$i});
			$add3 = 0; $add2 = 0; $add1 = 0;
			$add3 = $re2count{$re_loc}{"set3"} if exists $re2count{$re_loc}{"set3"}; $add3 = 2 if $add3 > 2;
			$add2 = $re2count{$re_loc}{"set2"} if exists $re2count{$re_loc}{"set2"}; $add2 = 2 if $add2 > 2;
			$add1 = $re2count{$re_loc}{"set1"} if exists $re2count{$re_loc}{"set1"}; $add1 = 2 if $add1 > 2;
			$count{$re_loc}{"set3"} += $add3;
			$count{$re_loc}{"set2"} += $add2;
			$count{$re_loc}{"set1"} += $add1;
			$all_3{$re_loc}{$i} += $add3;
			$all_2{$re_loc}{$i} += $add2;
			$all_1{$re_loc}{$i} += $add1;
		}
	}
	for $re_loc(sort keys %count){
		undef $c1; undef $c2; undef $c3; undef $fi1; undef $fi2; undef $fi3;
		$c1 = $count{$re_loc}{"set1"};
		$c2 = $count{$re_loc}{"set2"};
		$c3 = $count{$re_loc}{"set3"};
		$fi1 = $c1/($sample_no * 2);
		$fi2 = $c2/($sample_no * 2);
		$fi3 = $c3/($sample_no * 2);
		print O "$tp\t$re_loc\t$fi3\t$fi2\t$fi1\n";
	}
print "$tp done.\n";
}

for $re_loc(sort keys %all_1){
	$c1 = 0; $c2 = 0; $c3 = 0; $fi1 = 0; $fi2 = 0; $fi3 = 0;
	for $sample(keys %{$all_1{$re_loc}}){
		$add1 = 0; $add2 = 0; $add3 = 0;
		$add1 = $all_1{$re_loc}{$sample}; $add1 = 2 if $add1 > 2;
		$add2 = $all_2{$re_loc}{$sample}; $add2 = 2 if $add2 > 2;
		$add3 = $all_3{$re_loc}{$sample}; $add3 = 2 if $add3 > 2;
		$c1 += $add1; $c2 += $add2; $c3 += $add3;
	}
	$fi1 = $c1/($sample_no * 2);
	$fi2 = $c2/($sample_no * 2);
	$fi3 = $c3/($sample_no * 2);
	print O "all\t$re_loc\t$fi3\t$fi2\t$fi1\n";
}
close O || die;



