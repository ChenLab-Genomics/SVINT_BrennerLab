# calculate Frequency of Impact for the RE (FI), using the DGV SV set released on 6 March, 2019 (hg19; https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2_sv.sites.vcf.gz.tbi)
# selection criteria:
# (1) use all SVs: according to http://jmonlong.github.io/Hippocamplus/2019/03/31/gnomad-sv-first-look/, there are ~10,000 SVs that had significant overlaps with other SVs, but it was unclear whether they are real distinct ones called in different samples or duplicates due to technique issues; however, since 10,000 is only 2-3% of all gnomAD SVs, they are included anyways (as distinct ones)
# (2) only variants that are "PASS", or "MULTIALLELIC"; UNRESOLVED/PCRPLUS_ENRICHED/PREDICTED_GENOTYPING_ARTIFACT/VARIABLE_ACROSS_BATCHES were not included since they were undecided or likely technical effects
# SV types in this dataset
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=CPX,Description="Complex SV">
##ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INS:ME,Description="Mobile element insertion of unspecified ME class">
##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">
##ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">
##ALT=<ID=INS:ME:SVA,Description="SVA element insertion">
##ALT=<ID=INS:UNK,Description="Sequence insertion of unspecified origin">
##ALT=<ID=INV,Description="Inversion">
##MCNV
# (3) for insertions, the pos and end are used for their regions; BND and CTX are both counted as "translocations" (T), and CPX is counter as "others" (OTHER)
# (4) for insertions, count the affected region as the POS-1~POS (also, it's strange that in the VCF file, 768 INS had |END-POS| > 1Mb; it's unclear how the END were identified for insertions); their RE impacted always = 0.5 since insertions will introduce extra sequences but hard to evaluate their effects
# (5) for MCNV, the maf value is a sum of all types of <CN> except <CN=2>, which was the normal diplotype (presumably what most normal people have?)

use warnings;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require SV_SUB;

$sv_file = "data/gnomad/gnomad_v2_sv.sites.vcf";
$re_anno = "data/reg_anno/combined_regulatory_elements_for_SVint.bed"; #sorted

open F, $sv_file || die;
$tmp = <F>;
while(<F>){
	chomp; $line = $_;
	next if $line =~ /^#/;
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
	undef $id; undef $maf; undef $len; undef $chr; undef $stt; undef $end; undef $pos; undef $alt; undef $tp;
	($chr, $pos, $id, $ref, $alt, $qual, $fil, $info) = split /\t/, $line;
	next unless $fil =~ /PASS|MULTIALLELIC/;
	$stt = $pos - 1; 
	$chr = "chr$chr"; # original chr was 1,2,...
		@tmp = split ";", $info;
	undef $afs;
	for $t(@tmp){
		if($t =~ /END=(.*)/){
			$end = $1; next;
		}
		if($t =~ /SVTYPE=(.*)/){
			$tp = $1; next;
		}
		if($t =~ /^AF=(.*)/){
			$afs = $1; next;
		}
		if($t =~ /^SVLEN=(.*)/){
		$len = $1;
		}
	}
	if($tp =~ /MCNV/){
		$maf = 0; undef @alts; undef @afs; undef $i;
		@alts = split ",", $alt;
		@afs = split ",", $afs;
		for $i(0 .. $#afs){
			next if $alts[$i] =~ /CN\=2/;
			$af = $afs[$i];
			$maf = $maf + $af;
		}
	}else{
# print "$afs\n" if $afs =~ /,/;
		$maf = $afs;
	}
	if($stt > $end){
		$t = $stt; $stt = $end; $end = $t;
	}
# print "$maf\n" if $maf !~ /\d/;
	$var = "$chr\t$stt\t$end";
# define types 
	if($tp =~ /DEL/){
# DEL (deletion)
		if(!exists $to_use{"DEL"}{$var}){
			$to_use{"DEL"}{$var} = $maf;
		}else{ 
			$to_use{"DEL"}{$var} += $maf;
		}
		if(!exists $to_use{"all"}{$var}){
			$to_use{"all"}{$var} = $maf;
		}else{
			$to_use{"all"}{$var} += $maf;
		}
	}elsif($tp =~ /DUP/){
# DUP (duplication, tandem duplication)
		if(!exists $to_use{"DUP"}{$var}){
			$to_use{"DUP"}{$var} = $maf;
		}else{  
			$to_use{"DUP"}{$var} += $maf;
		}
		if(!exists $to_use{"all"}{$var}){
			$to_use{"all"}{$var} = $maf;
		}else{  
			$to_use{"all"}{$var} += $maf;
		}
	}elsif($tp =~ /MCNV/){
#MCNV (multi-allelic)
		if(!exists $to_use{"MCNV"}{$var}){
			$to_use{"MCNV"}{$var} = $maf;
		}else{  
			$to_use{"MCNV"}{$var} += $maf;
		}
		if(!exists $to_use{"all"}{$var}){
			$to_use{"all"}{$var} = $maf;
		}else{  
			$to_use{"all"}{$var} += $maf;
		}
	}elsif($tp =~ /INV/){
# INV
		if(!exists $to_use{"INV"}{$var}){
			$to_use{"INV"}{$var} = $maf;
		}else{  
			$to_use{"INV"}{$var} += $maf;
		}
		if(!exists $to_use{"all"}{$var}){
			$to_use{"all"}{$var} = $maf;
		}else{  
			$to_use{"all"}{$var} += $maf;
		}
	}elsif($tp =~ /INS/){
# INS (insertion, mobile element insertion, novel sequence insertion)
$posi = $pos - 1; $var2 = "$chr\t$posi\t$pos";
		if(!exists $to_use{"INS"}{$var2}){
			$to_use{"INS"}{$var2} = $maf;
		}else{  
			$to_use{"INS"}{$var2} += $maf;
		}
		if(!exists $to_use{"all"}{$var2}){
			$to_use{"all"}{$var2} = $maf;
		}else{  
			$to_use{"all"}{$var2} += $maf;
		}
	}elsif($tp =~ /BND|CTX/){
# BND, CTX: translocation
		if(!exists $to_use{"T"}{$var}){
			$to_use{"T"}{$var} = $maf;
		}else{  
			$to_use{"T"}{$var} += $maf;
		}
		if(!exists $to_use{"all"}{$var}){
			$to_use{"all"}{$var} = $maf;
		}else{  
			$to_use{"all"}{$var} += $maf;
		}
	}else{
# CPX
		if(!exists $to_use{"OTHER"}{$var}){
			$to_use{"OTHER"}{$var} = $maf;
		}else{  
			$to_use{"OTHER"}{$var} += $maf;
		}
		if(!exists $to_use{"all"}{$var}){
			$to_use{"all"}{$var} = $maf;
		}else{  
			$to_use{"all"}{$var} += $maf;
		}
	}
}
close F || die;

open T, ">tmp2.txt";
for $tp(sort keys %to_use){
	for $var(sort keys %{$to_use{$tp}}){
		print T "$tp\t$var\t$to_use{$tp}{$var}\n";
	}
}
close T;

&SV_SUB::write_bed2(\%to_use, "$sv_file\_temp");
$sv2re = SV_SUB::overlap("$sv_file\_temp.sorted.bed", $re_anno);
%sv2re = %{$sv2re};
unlink "$sv_file\_temp.sorted.bed";
unlink "$sv_file\_temp.bed";

# $..{$tp}{$sample_id}{$re}; for each sample, the same type of impact on the same RE will only be count once; for the same sample, different impact will all be taken into consideration
open O, ">$sv_file\_based_FI.txt" || die;
print O "SVTYPE\tCHR\tSTT\tEND\tFI(100% Affected)\tFI(50-100% Affected)\tFI(0-50% Affected)\n";
for $tp(keys %to_use){
	print "$tp\n";
	undef %re2fi;
	for $var(keys %{$to_use{$tp}}){
		$maf = $to_use{$tp}{$var};
		for $re_loc(keys %{$sv2re{$var}}){
			$ratio = $sv2re{$var}{$re_loc};
			undef $im;
			if($ratio >= 1){
				$im = "set3";
			}elsif($ratio < 0.5){
				$im = "set1";
			}else{
				$im = "set2";
			}
			if($tp eq "INS"){
				$im = "set2" if $im eq "set1";
			}
			if(exists $re2fi{$re_loc}{$im}){
				undef $c;
				$c = $re2fi{$re_loc}{$im};
				$re2fi{$re_loc}{$im} = $c + $maf;
			}else{
				$re2fi{$re_loc}{$im}=$maf;
			}
		}
	}
# 100% impact (set3) will also be counted for 0-0.5 (set1) and 0.5-1 (set2); 0.5-1 (set2) impact will also be counted for 0-0.5 (set1);
# for each case, the same re_loc with each impact should at most be counted twice (diplotype) 
	for $re_loc(keys %re2fi){
		$fi1 = 0; $fi2 = 0; $fi3 = 0;
		$fi1 = $re2fi{$re_loc}{"set1"} if exists $re2fi{$re_loc}{"set1"};
		$fi2 = $re2fi{$re_loc}{"set2"} if exists $re2fi{$re_loc}{"set2"};
		$fi3 = $re2fi{$re_loc}{"set3"} if exists $re2fi{$re_loc}{"set3"};
		$fi2 = $fi3 + $fi2;
		$fi1 = $fi1 + $fi2;
$fi1 = 1 if $fi1 > 1; $fi2 = 1 if $fi2 > 1; $fi3 = 1 if $fi3 > 1;
		print O "$tp\t$re_loc\t$fi3\t$fi2\t$fi1\n";
	}
}
close O || die;


