#! /usr/bin/perl
# the above line should be adjusted according to the actual path of perl, e.g. /usr/bin/perl
# currently only support SVs on chr1-22, chrX and chrY

package SV_SUB;
use warnings;

sub add {
	my ($x, $y) = @_;
	return $x + $y;
}

# read in family information
# ID parents/proband gender affected/unaffected sv-file-location
sub read_fam {
	my ($file) = @_;
	my ($proband, $father, $mother, %id2gender, %id2status, %id2svfile); 
	open I, $file;
	while(<I>){
		next if $_ =~ /^#/;
		chomp; ($id, $nm, $gender, $dis_status, $svfile) = split /\t/;
		$id2gender{$id} = $gender;
		$id2status{$id} = $dis_status;
		$id2svfile{$id} = $svfile;
		if($nm eq "proband"){
			$proband = $id; next;
		}elsif($nm eq "father"){
			$father = $id; next;
		}elsif($nm eq "mother"){
			$mother = $id; next;
		}else{
			push @relatives, $id;
		}
	}
	close I || die;
	return (\$proband, \$father, \$mother, \@relatives, \%id2gender, \%id2status, \%id2svfile);
}

# change TXT input into VCF input
# TXT input format: chr, stt, end, id, type
sub write_txt2vcf{
	my ($file) = @_;
	open I, $file;
	open O, ">$file".".fake.vcf";
	print O "##fileformat=VCFv4.1\n##this file is converted from $file into VCF format\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
	while(<I>){
		chomp; $_ =~ s/\r$//; # in some OS txt files would have \r at the end of each line...\r need to be removed
		next if $_ !~ /^chr\d|^chrX|^chrY/i; 
		undef $tmp; $tmp = $_;
		undef $chr; undef $stt; undef $end; undef $id; undef $tp; undef $pos; undef $ref; undef $alt; undef $fil; undef $info; undef $len; undef $format; undef $gt;
		($chr, $stt, $end, $id, $tp) = split /\t/, $tmp;
		$chr =~ s/Chr|CHR/chr/ if $chr =~ /Chr|CHR/;
		$pos = $stt; # ignore the 0-based or 1-based issue: would not differ much for an SV
			$len = $end - $stt;
		$ref = "."; $alt = $tp; $qual = "."; $fil = "PASS";
		$tp =~ s/\<// if $tp =~ /\</; $tp =~ s/\>// if $tp =~ /\>/;
		$info = "SVTYPE=$tp;END=$end;SVLEN=$len"; $format = "GT"; $gt = "0/1";
		print O "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$fil\t$info\t$format\t$gt\n";

	}
	close I || die;
	close O || die;
	return ();
}

# read in SV files in VCF format (using 10X longranger call format as an example)
# output: same chr SV; cross-chr SV
# use only paired BND; use "0/1" if two BND disagree on GT
sub read_sv {
	my ($file) = @_;
	my (%same, %diff, %mates, %bnd_lines, %lines);
	open I, $file;
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BC00101_longranger_noloupe_GATK
	my $count=0;	
	while(<I>){
		next if $_ =~ /^#/;
		chomp; $line=$_; ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $tgt) = split /\t/;
		$count++;
		$bnd_lines{$id}{$line} = 1; $id2count{$id}=$count; $lines{$count}{$line}=1;
		$chr=~ s/Chr/chr/; $chr = "chr$chr" if $chr !~ /chr/;
		print "warning for SV call: format not contain \"GT\"!\n" if $format !~ "GT";
		@f = split ":", $format; @tgt = split ":", $tgt;
		for $k(0 .. $#f){
			if($f[$k] eq "GT"){
				$gt = $tgt[$k]; last;
			}
		}
		if($alt =~ /^<(.*)>$/){
			$tp = $1;
			$end = $pos + 1;
			if($info =~ /END=(\d+)/){
				$end = $1;
			}elsif($info =~ /SVLEN=(\d+);/){
				$end = $pos + $1;
			}
			@loc = sort{$a<=>$b}($pos, $end); $var = "$chr\t$loc[0]\t$loc[1]";
			$same{$var}{$tp} = $gt; $lines{$count}{$line}=$var;
		}else{
			undef $m_id; $e=0;
			@in = split ";", $info;
			for $in(@in){
				if($in =~ /MATEID=(.*)/){
					$m_id = $1; $e=1;
					last;
				}}
			if($e == 1){
				@ids = sort{$a cmp $b}($id, $m_id);
				$mates{$ids[0]}{$ids[1]}=1;
			}
		}
	}
	close I || die;
	for $id1(keys %mates){
		next if !exists $id2count{$id1};
		@line1 = keys %{$bnd_lines{$id1}}; @t1 = split /\t/, $line1[0];$count1=$id2count{$id1};
		for $id2(keys %{$mates{$id1}}){
			next if !exists $id2count{$id2};
			@line2 = keys %{$bnd_lines{$id2}}; @t2 = split /\t/, $line2[0]; $count2=$id2count{$id2};
			$chr1 = $t1[0]; $pos1 = $t1[1]; $chr2 = $t2[0]; $pos2 = $t2[1]; $gt1=$t1[$#t1]; $gt2=$t2[$#t2];
			$chr1 =~ s/Chr/chr/; $chr1 = "chr$chr1" if $chr1 !~ /chr/; $chr2 =~ s/Chr/chr/; $chr2 = "chr$chr2" if $chr2 !~ /chr/;
			undef $tp; undef $gt;
			if($gt1 ne $gt2){
				$gt="0/1";
			}else{
				$gt=$gt1;
			}
			@in = split ";", $t1[7];
			for $in(@in){
				if($in =~ /SVTYPE2=(.*)/){
					$tp = $1;
				}elsif($in =~ /SVTYPE=(.*)/){
					$tp = $1;
				}
			}
			if($chr1 eq $chr2){
				@loc = sort{$a<=>$b}($pos1, $pos2); $var = "$chr1\t$loc[0]\t$loc[1]";
				$same{$var}{$tp}=$gt; 
				$lines{$count1}{$line1[0]}=$var; $lines{$count2}{$line2[0]}=$var;
			}else{
				$var1 = "$chr1\t$pos1"; $var2 = "$chr2\t$pos2";
				@loc = sort{$a cmp $b}($var1, $var2); $var = "$loc[0]\t$loc[1]";
				$diff{$var}{$tp}=$gt;
				$lines{$count1}{$line1[0]} = $var; $lines2{$count2}{$line2[0]} = $var;
			}
		}}
	return(\%same, \%diff, \%lines);
}


#
sub write_bed{ # based on the result from read_in()
	my ($list, $bed) = @_;
	my %list = %{$list};
	open O, ">$bed.bed";
	for $var(keys %list){
		($chr, $stt, $end) = split /\t/, $var;
		print O "$chr\t$stt\t$end\n";
	}
	close O || die;
	system "sort -k1,1 -k2,2n $bed.bed > $bed.sorted.bed";
	return()
}

#
sub write_bed2{ # based on the result from read_in()
	my ($list, $bed) = @_;
	my %list = %{$list};
	open O, ">$bed.bed";
	for $tp(keys %list){
		for $var(keys %{$list{$tp}}){
			($chr, $stt, $end) = split /\t/, $var;
			print O "$chr\t$stt\t$end\n";
		}        
	}
	close O || die;
	system "sort -k1,1 -k2,2n $bed.bed > $bed.sorted.bed";
	return()
}

#
sub overlap{ # overlap("query.bed", "feature.bed"); need to be sorted
	my ($query, $feature) = @_;
	my %over;
	open I, "bedtools intersect -wa -wb -a $query -b $feature -sorted |";
	while(<I>){
		chomp;
		($chr1, $stt1, $end1, $chr2, $stt2, $end2) = split /\t/;
		@tmp = sort{$a <=> $b}($stt1, $end1, $stt2, $end2);
		$over{"$chr1\t$stt1\t$end1"}{"$chr2\t$stt2\t$end2"} = ($tmp[2] - $tmp[1])/($end2-$stt2);
	}
	close I;
	return (\%over);
} 

#
sub overlap2{ # overlap("query.bed", "feature.bed"); need to be sorted
	my ($query, $feature) = @_;
	my %over;
	open I, "bedtools intersect -wo -a $query -b $feature -sorted |";
	while(<I>){
		chomp;
		($chr1, $stt1, $end1, $chr2, $stt2, $end2, $gn) = split /\t/;
		@tmp = sort{$a <=> $b}($stt1, $end1, $stt2, $end2);
		$over{"$chr1\t$stt1\t$end1"}{$gn} = 1;
	}
	close I;
	return (\%over);
}

#
sub anno_tad1{ # anno_tad1("query.bed", "boundary_regions.bed", "data/tad/cell_bd_nearby_genes.txt"); need to be sorted
	my ($query, $feature, $feature2) = @_;
	my %over; my %ci2gn;
	open F, $feature2;
	while(<F>){
		chomp; ($celli, $gns) = split /\t/;
		if($gns){ #some boundaries don't have genes nearby
		$ci2gn{$celli} = $gns;
		}
	}
	close F;
	open I, "bedtools intersect -wo -a $query -b $feature -sorted |";
	while(<I>){
		chomp; undef $cell;
		($chr1, $stt1, $end1, $chr2, $stt2, $end2, $celli, $t, $tt, $bp) = split /\t/;
		print "TAD boundary region name error: $celli\n" if $celli !~ /\:/;
		($cell, $i) = split ":", $celli;
		$add = "-";
		$add = $ci2gn{$celli} if exists $ci2gn{$celli};
		if(!exists $over{"$chr1\t$stt1\t$end1"}{$cell}){
			$over{"$chr1\t$stt1\t$end1"}{$cell} = $add;
		}else{
			$over{"$chr1\t$stt1\t$end1"}{$cell} = $over{"$chr1\t$stt1\t$end1"}{$cell}.",".$add;		
		}
	}
	close I;
	for $var(keys %over){
		for $cell(keys %{$over{$var}}){
			undef %gn;
			for $g(split ",", $over{$var}{$cell}){
				$gn{$g} = 1 if $g ne "-"; #some boundaries don't have genes nearby
			}
			if(scalar(keys %gn) > 0){
				$over{$var}{$cell} = join ",", sort keys %gn;
			}else{
				$over{$var}{$cell} = "-";
			}
		}
	}
	return (\%over);
}

# annotate SV with RE, targets, and impact group of RE
sub rg_anno{
	my ($sv2re, $re_anno, $rg_file, $overlap_gene) = @_;
	my %sv2re = %{$sv2re}; my %overlap_gene = %{$overlap_gene};
	my %loc2re; my %sv2target; my %loc_gn;
	open F, $re_anno;
	while(<F>){
		($chr, $stt, $end, $re) = split /\t/;
		$loc2re{"$chr\t$stt\t$end"}=$re;
	}
	close F || die;
	open I, $rg_file || die;
	while(<I>){
		chomp; ($re, $gn, $ms) = split /\t/; # regulatory element; gene; methods
			$re2gn{$re}{$gn}=$ms;
	}
	close I || die;
	for $var(keys %sv2re){
		my @set1; my @set2; my @set3;
		if(scalar keys %{$overlap_gene{$var}}>=1){ # genic noncoding region SVs
			for $gn(keys %{$overlap_gene{$var}}){
				for $re_loc(keys %{$sv2re{$var}}){
					$ratio = $sv2re{$var}{$re_loc};
					$re = $loc2re{$re_loc}; $loc = join ":", (split "\t", $re_loc);
					if($ratio < 0.5){
						push @set1, $loc."--".$gn."(Overlapped)";
						$loc_gn{$loc}{$gn}=1;
					}elsif($ratio >= 1){
						push @set3, $loc."--".$gn."(Overlapped)";
						$loc_gn{$loc}{$gn}=1;						
					}else{
						push @set2, $loc."--".$gn."(Overlapped)";
						$loc_gn{$loc}{$gn}=1;				
					}

				}
			}
		}else{ # intergenic SVs
			for $re_loc(keys %{$sv2re{$var}}){
				$ratio = $sv2re{$var}{$re_loc};
				$re = $loc2re{$re_loc}; $loc = join ":", (split "\t", $re_loc);
				if($ratio < 0.5){
					for $gn(keys %{$re2gn{$re}}){
						push @set1, $loc."--".$gn."($re2gn{$re}{$gn})";
						$loc_gn{$loc}{$gn}=1;
					}
				}elsif($ratio >= 1){
					for $gn(keys %{$re2gn{$re}}){
						push @set3, $loc."--".$gn."($re2gn{$re}{$gn})";
						$loc_gn{$loc}{$gn}=1;				
					}			
				}else{
					for $gn(keys %{$re2gn{$re}}){
						push @set2, $loc."--".$gn."($re2gn{$re}{$gn})";
						$loc_gn{$loc}{$gn}=1;				
					}
				}
			}
		}
		my $j3 = join ";", @set3; $j3 = "-" if $j3 eq "";
		my $j2 = join ";", @set2; $j2 = "-" if $j2 eq "";
		my $j1 = join ";", @set1; $j1 = "-" if $j1 eq "";
		$sv2target{$var}="$j3|$j2|$j1";
	}
	return (\%sv2target, \%loc_gn);
}

#
sub write_bed3{ # based on the result from rg_anno()
	my ($loc_gn, $gn_anno, $bed) = @_;
	my %loc_gn = %{$loc_gn};
	my %gn2tss;
	@chr=("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY", "chrM", "chrMT"); map{$chr{$_}=1}@chr;
	open G, $gn_anno;
	while(<G>){
		chomp; @t = split /\t/, $_;
		next if !exists $chr{$t[0]};
		if($t[5] eq "+"){
			$gn2tss{$t[3]} = $t[1];
		}else{
			$gn2tss{$t[3]} = $t[2];
		}
		$gn2chr{$t[3]}=$t[0];
	}
	close G || die;
	open O, ">$bed.bed";
	for $loc(keys %loc_gn){
		for $gn(keys %{$loc_gn{$loc}}){
			($chr, $stt, $end) = split /\:/, $loc;
			if (!exists $gn2tss{$gn}){
				$tss = $stt; # some genes are missing from hg38 gtf; for such genes, there are no way to find out whether the RE-gene interaction is across TAD boundaries; only the RE region were tested for TAD boundaries
			}else{ 
				$tss = $gn2tss{$gn};
			}
			if(exists $gn2chr{$gn}){
				$g_chr = $gn2chr{$gn};
			}else{
				$g_chr = $chr;  # some genes are missing from hg38 gtf; for such genes, there are no way to find out whether the RE-gene interaction is across TAD boundaries; only the RE region were tested for TAD boundaries
			}
			if($chr eq $g_chr){
				@tt = ($tss, $stt, $end);
				@tt = sort{$a<=>$b}@tt;
				print O "$chr\t$tt[0]\t$tt[2]\t$loc--$gn\n";
			}else{
				print O "$chr\t0\t3000000000\t$loc--$gn\n";
			}
		}
	}
	close O || die;
	system "sort -k1,1 -k2,2n $bed.bed > $bed.sorted.bed";
	return()
}

#
sub anno_tad2{ # anno_tad1("from-write_bed3.bed", "boundary_regions.bed"); need to be sorted
	my ($outfile, $query, $feature) = @_;
	my %over; my %over2;
	open T, $query; open OT, ">$query.tmp";
	while(<T>){
		chomp; $t=$_;
		($chr, $stt, $end, $lg) = split /\t/;
		if($end-$stt>10000000){
			$over{$lg} = "Distance>10Mb";
		}else{
			print OT $t."\n";
		}
	}
	close T || die; close OT || die;
	system "sort -k1,1 -k2,2n $query.tmp > $query.tmp.sorted";
	open I, "bedtools intersect -wao -a $query.tmp.sorted -b $feature -sorted |";
	while(<I>){
		chomp; undef $cell;
		($chr1, $stt1, $end1, $lg, $chr2, $stt2, $end2, $celli, $t, $tt, $bp) = split /\t/;
		if($bp == 0){
			$over{$lg} = 0;
		}else{
			print "TAD boundary region name error: $celli\n" if $celli !~ /\:/;
			($cell, $i) = split ":", $celli;
			$over2{$lg}{$cell} = 1;
		}
	}
	close I;
	open O, ">$outfile";
	for $lg(keys %over){
		if($over{$lg} eq "0"){
			print O "$lg\tNone\n"; next;
		}
		print O "$lg\t$over{$lg}\n";
	}
	for $lg(keys %over2){
		@cs = keys %{$over2{$lg}};
		$cs = join ",", @cs;
		print O "$lg\t$cs\n";
		$over{$lg} = scalar @cs;
	}
	close O || die;
	return (\%over);
}

# annotate FI for RE based on 1000Genome SV calls
sub read_fi{
	my ($file, $sv2re) = @_;
	my %fi;
	my %sv2re = %{$sv2re};
	my %see_re;
	for $sv(keys %sv2re){
	for $re_loc(keys %{$sv2re{$sv}}){
		$see_re{$re_loc} = 1;
		}
	}
	open I, $file || die; 
	$tmp=<I>;
	while(<I>){
		chomp;
		($tp, $chr, $stt, $end, $fi3, $fi2, $fi1)=split /\t/;
		$re_loc = "$chr:$stt:$end"; $re_loc2 = "$chr\t$stt\t$end";
		next if(!exists $see_re{$re_loc2});
		$fi{$tp}{$re_loc}{"set3"} = sprintf("%.3e", $fi3);
		$fi{$tp}{$re_loc}{"set2"} = sprintf("%.3e", $fi2);
		$fi{$tp}{$re_loc}{"set1"} = sprintf("%.3e", $fi1);
	}
	close I || die;
	return (\%fi);
}

sub anno_fi_tad{
	my ($sv2target, $to_use, $fi1, $fi2, $lg2tad_bd) = @_;
	%fi1 = %{$fi1}; %fi2 = %{$fi2}; %lg2tad_bd = %{$lg2tad_bd}; %sv2target = %{$sv2target}; %to_use = %{$to_use}; 
	my %sv2fi_tad;
	for $var(keys %sv2target){
		next if(!exists $to_use{$var});
		for $tp0(keys %{$to_use{$var}}){
			if($tp0 =~ /DUP/){
				$tp = "DUP";
			}elsif($tp0 =~ /DEL/){
				$tp = "DEL";
			}elsif($tp0 =~ /CNV/){
				$tp = "CNV";
			}elsif($tp0 =~ /INV/){
				$tp = "INV";
			}elsif($tp0 =~ /INS/){
				$tp = "INS";
			}else{
				$tp = "all"; # if cannot match the SV type, just use the FI by all SVs impacting the RE
			}
			undef @ts3; undef @ts2; undef @ts1; my @ts = ("-","-","-");
			($ts[0], $ts[1], $ts[2]) = split /\|/, $sv2target{$var};
			undef @set3; undef @set2; undef @set1;
			if($ts[0] ne "-"){
				@ts3 = split ";", $ts[0];
				for $tmp(@ts3){
					undef $lg; undef $ms; undef $re_loc; undef $f; undef @cs; undef $tad;
					$tmp =~ /(.*)\((.*)\)/; $lg = $1; $ms = $2;
					if($lg2tad_bd{$lg} =~ /Distance/){
						$tad = $lg2tad_bd{$lg};
					}else{
						$tad = "SameTAD=".(22-$lg2tad_bd{$lg})."/22";					
					}
					($re_loc) = split "--", $lg; # $lg is "re_loc--gn"
					undef $f1; undef $f2;
						if(exists $fi1{$tp}{$re_loc}{"set3"}){
							$f1= "1000G_FI=".$fi1{$tp}{$re_loc}{"set3"}."[$tp]";
# if no FI for the specific SV type, try using the FI by all SVs
						}elsif(exists $fi1{"all"}{$re_loc}{"set3"}){
							$f1= "1000G_FI=".$fi1{"all"}{$re_loc}{"set3"}."[all]";
						}else{
							$f1= "1000G_FI=0[all]";
						}
						if(exists $fi2{$tp}{$re_loc}{"set3"}){
                                                        $f2= "gnomAD_FI=".$fi2{$tp}{$re_loc}{"set3"}."[$tp]";
# if no FI for the specific SV type, try using the FI by all SVs
                                                }elsif(exists $fi2{"all"}{$re_loc}{"set3"}){
                                                        $f2= "gnomAD_FI=".$fi2{"all"}{$re_loc}{"set3"}."[all]";
                                                }else{
                                                        $f2= "gnomAD_FI=0[all]";
                                                }
					push @set3, $lg."($ms!$f1,$f2!$tad)";
				}
			}
			if($ts[1] ne "-"){
				@ts2 = split ";", $ts[1];
				for $tmp(@ts2){
					undef $lg; undef $ms; undef $re_loc; undef $f; undef @cs; undef $tad;
					$tmp =~ /(.*)\((.*)\)/; $lg = $1; $ms = $2;
					if($lg2tad_bd{$lg} =~ /Distance/){
						$tad = $lg2tad_bd{$lg};
					}else{
						$tad = "SameTAD=".(22-$lg2tad_bd{$lg})."/22";           
					}					
					($re_loc) = split "--", $lg; # $lg is "re_loc--gn"
					undef $f1; undef $f2;
					if(exists $fi1{$tp}{$re_loc}{"set3"}){
						$f1= "1000G_FI=".$fi1{$tp}{$re_loc}{"set3"}."[$tp]";
					}elsif(exists $fi1{"all"}{$re_loc}{"set3"}){
						$f1= "1000G_FI=".$fi1{"all"}{$re_loc}{"set3"}."[all]";
					}else{
						$f1= "1000G_FI=0[all]";
					}
					if(exists $fi2{$tp}{$re_loc}{"set3"}){
						$f2= "gnomAD_FI=".$fi2{$tp}{$re_loc}{"set3"}."[$tp]";
					}elsif(exists $fi2{"all"}{$re_loc}{"set3"}){
						$f2= "gnomAD_FI=".$fi2{"all"}{$re_loc}{"set3"}."[all]";
					}else{  
						$f2= "gnomAD_FI=0[all]";
					}
					push @set2, $lg."($ms!$f1,$f2!$tad)";
				}
			}
			if($ts[2] ne "-"){
				@ts1 = split ";", $ts[2];
				for $tmp(@ts1){
					undef $lg; undef $ms; undef $re_loc; undef $f; undef @cs; undef $tad;
					$tmp =~ /(.*)\((.*)\)/; $lg = $1; $ms = $2;
#print "$lg\:$lg2tad_bd{$lg}\n";					
					if($lg2tad_bd{$lg} =~ /Distance/){
						$tad = $lg2tad_bd{$lg};
					}else{
						$tad = "SameTAD=".(22-$lg2tad_bd{$lg})."/22";           
					}					
					($re_loc) = split "--", $lg; # $lg is "re_loc--gn"
					undef $f1; undef $f2;
					if(exists $fi1{$tp}{$re_loc}{"set3"}){
						$f1= "1000G_FI=".$fi1{$tp}{$re_loc}{"set3"}."[$tp]";
					}elsif(exists $fi1{"all"}{$re_loc}{"set3"}){
						$f1= "1000G_FI=".$fi1{"all"}{$re_loc}{"set3"}."[all]";
					}else{
						$f1= "1000G_FI=0[all]";
					}
					if(exists $fi2{$tp}{$re_loc}{"set3"}){
						$f2= "gnomAD_FI=".$fi2{$tp}{$re_loc}{"set3"}."[$tp]";
					}elsif(exists $fi2{"all"}{$re_loc}{"set3"}){
						$f2= "gnomAD_FI=".$fi2{"all"}{$re_loc}{"set3"}."[all]";
					}else{  
						$f2= "gnomAD_FI=0[all]";
					}
					push @set1, $lg."($ms!$f1,$f2!$tad)";
				}
			}
			my $j3 = join ";", @set3; $j3 = "-" if $j3 eq "";
			my $j2 = join ";", @set2; $j2 = "-" if $j2 eq "";
			my $j1 = join ";", @set1; $j1 = "-" if $j1 eq "";
			$sv2fi_tad{$var}{$tp0} = "$j3|$j2|$j1"; 
		}
	}
	return (\%sv2fi_tad);
}

#
sub write_svint_anno{
	my ($out_file, $sv_file, $lines_a, $same_a, $to_use, $overlap_gene, $overlap_exon, $sv2tad_bd, $sv2fi_tad) = @_;
	%sv2fi_tad = %{$sv2fi_tad}; %sv2tad_bd=%{$sv2tad_bd}; %lines_a = %{$lines_a}; %same_a = %{$same_a}; %to_use = %{$to_use}; %overlap_gene = %{$overlap_gene}; %overlap_exon = %{$overlap_exon};
	open O, ">$out_file" || die;
	print O "##SVint annotated target genes, TAD boundary impact, and FI (Frequency of Impact for Regulatory elements)  for noncoding SVs\n";
	open I, $sv_file || die;
	while(<I>){
		chomp; $tmp = $_; last if $tmp !~ /^#/;
		if($tmp =~ /^##/){
			print O $tmp."\n";
		}elsif($tmp =~ /^#CHR/){
			print O "##SVINT Format: Celllines-where-SV-Impact-TAD-boundary-regions|regulatory regions affected 100%|regulatory regions affected 50-100%|regulatory regions affected 0-50%\n##SVINT Format details in each of the three affected sections: RegulatoryRegion--targetGene(prediction methods!Frequency-of-Impact!No of celllines in which RegulatoryRegion--targetGene are in same TADs)\n"; 
			print O "$tmp\tOverlapWithGenes\tSVINT\n";
		}
	}
	close I || die;
	for $id(sort{$a<=>$b}(keys %lines_a)){
		for $line(keys %{$lines_a{$id}}){
			undef @svint; undef @tt; undef $tp0;
			$var = $lines_a{$id}{$line};
			if(!exists $same_a{$var}){
				print O "$line\tCross different chromosomes\t-|-|-|-\n"; 
				next;
			}
			if(scalar keys %{$to_use{$var}}<1){ # exon overlapped SVs
				undef $ogs; $ogs = join ",", (keys %{$overlap_gene{$var}});
				if(scalar keys %{$overlap_exon{$var}}>=1){
					print O "$line\tOverlap with exonic regions of genes: $ogs\t-|-|-|-\n"; 
					next;
				}
			}
			@tt = split /\t/, $line; $alt = $tt[4];
			if($alt =~ /<(.*)>/){
				$tp0 = $1;
			}else{
				undef @in;
				@in = split ";", $tt[7];
				for $in(@in){
					if($in =~ /SVTYPE2=(.*)/){
						$tp0 = $1;
					}elsif($in =~ /SVTYPE=(.*)/){
						$tp0 = $1;
					}
				}
			}
			undef @cs; @cs = keys %{$sv2tad_bd{$var}};
			if((scalar @cs)>0){
				undef @cg;
				for $c(@cs){
					push @cg, $c."(".$sv2tad_bd{$var}{$c}.")"; # cell1(affected genes 1);cell2(affected genes 2)...
				}
				push @svint, (join ";", sort @cg);
			}else{
				push @svint, "-";
			}
			if(scalar keys %{$overlap_gene{$var}}>=1){ # genic noncoding region SVs
				undef $ogs; $ogs = join ",", (keys %{$overlap_gene{$var}});
				if(!exists $sv2fi_tad{$var}{$tp0}){ # no regulatory elements overlapped
					print O "$line\tnoncoding regions of genes: $ogs\t$svint[0]|-|-|-\n"; 
					next;
				}else{
					push @svint, $sv2fi_tad{$var}{$tp0};
					print O "$line\tnoncoding regions of genes: $ogs\t".(join "|", @svint)."\n";
					next;
				}
			}else{ # intergenic SVs
				if(!exists $sv2fi_tad{$var}{$tp0}){ # no regulatory elements overlapped
					push @svint, "-|-|-";
				}else{
					push @svint, $sv2fi_tad{$var}{$tp0};
				}
				print O "$line\tNone\t".(join "|", @svint)."\n";
			}
		}
	}
	close O || die;
}

#
sub write_svint_anno_incluEx{
	my ($out_file, $sv_file, $lines_a, $same_a, $to_use, $overlap_gene, $overlap_exon, $gn2ex, $sv2tad_bd, $sv2fi_tad) = @_;
	%sv2fi_tad = %{$sv2fi_tad}; %sv2tad_bd=%{$sv2tad_bd}; %lines_a = %{$lines_a}; %same_a = %{$same_a}; %to_use = %{$to_use}; %overlap_gene = %{$overlap_gene}; %overlap_exon = %{$overlap_exon}; %gn2ex = %{$gn2ex};
	open O, ">$out_file" || die;
	print O "##SVint annotated target genes, TAD boundary impact, and FI (Frequency of Impact for Regulatory elements)  for noncoding SVs\n";
	open I, $sv_file || die;
	while(<I>){
		chomp; $tmp = $_; last if $tmp !~ /^#/;
		if($tmp =~ /^##/){
			print O $tmp."\n";
		}elsif($tmp =~ /^#CHR/){
			print O "##SVINT Format: Celllines-where-SV-Impact-TAD-boundary-regions|regulatory regions affected 100%|regulatory regions affected 50-100%|regulatory regions affected 0-50%\n##SVINT Format details in each of the three affected sections: RegulatoryRegion--targetGene(prediction methods!Frequency-of-Impact!No of celllines in which RegulatoryRegion--targetGene are in same TADs)\n";
			print O "$tmp\tOverlapWithGenes\tSVINT\n";
		}
	}
	close I || die;
	for $id(sort{$a<=>$b}(keys %lines_a)){
		for $line(keys %{$lines_a{$id}}){
			undef @svint; undef @tt; undef $tp0;
			$var = $lines_a{$id}{$line};
			if(!exists $same_a{$var}){
				print O "$line\tCross different chromosomes\t-|-|-|-\n";
				next;
			}
#if(scalar keys %{$to_use{$var}}<1){ # exon overlapped SVs
#       undef $ogs; $ogs = join ",", (keys %{$overlap_gene{$var}});
#      if(scalar keys %{$overlap_exon{$var}}>=1){
#             print O "$line\,Overlap with exonic regions of genes: $ogs\t-|-|-|-\n";
#            next;
#   }
#}
			@tt = split /\t/, $line; $alt = $tt[4];
			if($alt =~ /<(.*)>/){
				$tp0 = $1;
			}else{
				undef @in;
				@in = split ";", $tt[7];
				for $in(@in){
					if($in =~ /SVTYPE2=(.*)/){
						$tp0 = $1;
					}elsif($in =~ /SVTYPE=(.*)/){
						$tp0 = $1;
					}
				}
			}
			undef @cs; @cs = keys %{$sv2tad_bd{$var}};
			if((scalar @cs)>0){
				undef @cg;
				for $c(@cs){
					push @cg, $c."(".$sv2tad_bd{$var}{$c}.")"; # cell1(affected genes 1);cell2(affected genes 2)...
				}
				push @svint, (join ";", sort @cg);
			}else{
				push @svint, "-";
			}
			if(scalar keys %{$overlap_gene{$var}}>=1){ # genic region SVs
				undef @ogs1; undef @ogs2; undef $ogs1; undef $ogs2;
				for $g(keys %{$overlap_gene{$var}}){
					$tmp = 0;
					for $e(keys %{$gn2ex{$g}}){
#print "$e\n";
						if(exists $overlap_exon{$var}{$e}){
							$tmp ++;
						}
					}
					if($tmp > 0){
						push @ogs1, $g; # exonic
					}else{
						push @ogs2, $g; # intronic
					}
				}
				$ogs1 = join ",", @ogs1; $ogs1 = "." if !ogs1;
				$ogs2 = join ",", @ogs2; $ogs2 = "." if !ogs2;
				if(!exists $sv2fi_tad{$var}{$tp0}){ # no regulatory elements overlapped
					print O "$line\tOverlap with exonic regions of genes: $ogs1; Overlap with only noncoding regions of genes: $ogs2\t$svint[0]|-|-|-\n";
					next;
				}else{
					push @svint, $sv2fi_tad{$var}{$tp0};
					print O "$line\tOverlap with exonic regions of genes: $ogs1; Overlap with only noncoding regions of genes: $ogs2\t".(join "|", @svint)."\n";
					next;
				}
			}else{ # intergenic SVs
				if(!exists $sv2fi_tad{$var}{$tp0}){ # no regulatory elements overlapped
					push @svint, "-|-|-";
				}else{
					push @svint, $sv2fi_tad{$var}{$tp0};
				}
				print O "$line\tNone\t".(join "|", @svint)."\n";
			}
		}
	}
	close O || die;
}


# write out TXT format results in several levels: FI=0 & 100% impact & >=2 methods; FI<=0.01 & 50-100% impact & >=2 methods; FI<=0.05 & 50-100% impact & >=1 methods;
sub write_txt_anno{
	my ($sv_anno_file, $fi_cut, $out_file1, $out_file2, $out_file3, $out_file4) =  @_;
# $out_file1: FI<=cut-off & 100% impact & >=2 methods
# $out_file2: FI<=cut-off & 100% impact & >=1 method
# $out_file2: FI<=cut_off & 50-100% impact & >=2 methods
# $out_file4: FI<=cut_off & 50-100% impact & >=1 method
	open I, $sv_anno_file;
	open O1, ">$out_file1"; 
	print O1 "#subset results from $sv_anno_file: for regulatory target genes, use max FI of RE (among different databases)<=$fi_cut & regulatory element region 100% impacted & supported by >=2 methods (if two methods supported two different regulatory element but the same target gene, this gene is also included)\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tOverlap_With_Gene_Exons\tGenic_Noncoding\tRegulatory_Target_Genes\tTAD_Boundary_Cellines\n";
	open O2, ">$out_file2";
	print O2 "#subset results from $sv_anno_file: for regulatory target genes, use max FI of RE (among different databases)<=$fi_cut & regulatory element region 100% impacted & supported by >=1 method\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tOverlap_With_Gene_Exons\tGenic_Noncoding\tRegulatory_Target_Genes\tTAD_Boundary_Cellines\n";
	open O3, ">$out_file3";
	print O3 "#subset results from $sv_anno_file: for regulatory target genes, use max FI of RE (among different databases)<=$fi_cut & regulatory element region >=50% impacted & supported by >=2 methods (if two methods supported two different regulatory element but the same target gene, this gene is also included)\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tOverlap_With_Gene_Exons\tGenic_Noncoding\tRegulatory_Target_Genes\tTAD_Boundary_Cellines\n";
	open O4, ">$out_file4";
	print O4 "#subset results from $sv_anno_file: for regulatory target genes, use max FI of RE (among different databases)<=$fi_cut & regulatory element region >=50% impacted & supported by >=1 method\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tOverlap_With_Gene_Exons\tGenic_Noncoding\tRegulatory_Target_Genes\tTAD_Boundary_Cellines\n";
	while(<I>){
		next if $_ =~ /^#/; chomp;
		undef $chr; undef $pos; undef $id; undef $ref; undef $alt; undef $qual; undef $filter; undef $info; undef $format; undef $gt; undef $over; undef $anno;
		($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $gt, $over, $anno) = split /\t/;
		$exon = "."; $noncoding = ".";
		for $o(split ";", $over){
			if($o =~ /exonic regions.*: (.*)/){
				$exon = $1; 
			}
			if($o =~ /noncoding regions.*: (.*)/){
				$noncoding = $1;
			}
		}
		undef $tad; undef $tar4; undef $tar3; undef $tar2; undef $tar1; undef %tar4; undef %tar3; undef %tar2; undef %tar1;
		if($anno =~ /\|/){
			($tad, $anno100, $anno50, $anno0) = split /\|/, $anno;
			$tad = "." if $tad eq "-";
			if($anno100 ne "-"){
				for $t(split ";", $anno100){
					undef $g; undef $fs; undef $ms;
					undef $f1; undef $f2; undef $fi1; undef $fi2;
					$t =~ /--(.*)\((.*)!(.*FI.*)!.*\)/;
					$g = $1; $ms = $2; $fs = $3;
					($f1, $f2) = split ",", $fs;
					if($f1 =~ /1000G_FI=(.*)\[/){
						$fi1 = $1;
					}
					if($f2 =~ /gnomAD_FI=(.*)\[/){
                                                $fi2 = $1;
                                        }
					next if $fi1 > $fi_cut;
					next if $fi2 > $fi_cut;
					$tar2{$g}=1;
					$tar4{$g}=1;
					for $m(split ",", $ms){
						$tar1{$g}{$m}=1;
						$tar3{$g}{$m}=1;
					}
				}
			}
			if($anno50 ne "-"){
				for $t(split ";", $anno50){
					undef $g; undef $fs; undef $ms;
                                        undef $f1; undef $f2; undef $fi1; undef $fi2;
                                        $t =~ /--(.*)\((.*)!(.*FI.*)!.*\)/;
                                        $g = $1; $ms = $2; $fs = $3;
                                        ($f1, $f2) = split ",", $fs;
                                        if($f1 =~ /1000G_FI=(.*)\[/){
                                                $fi1 = $1;
                                        }
                                        if($f2 =~ /gnomAD_FI=(.*)\[/){
                                                $fi2 = $1;
                                        }
                                        next if $fi1 > $fi_cut;
                                        next if $fi2 > $fi_cut;
					$tar4{$g}=1;
					for $m(split ",", $ms){
						$tar3{$g}{$m}=1;
					}
				}
			}
			for $g1(sort keys %tar1){
				next if scalar(keys %{$tar1{$g1}}) < 2;
				if(!$tar1){
					$tar1 = $g1;
				}else{
					$tar1 = "$tar1,$g1";
				}
			}
			for $g3(sort keys %tar3){
				next if scalar(keys %{$tar3{$g3}}) < 2;
				if(!$tar3){
					$tar3 = $g3;
				}else{
					$tar3 = "$tar3,$g3";
				}
			}
			$tar2 = join ",", sort keys %tar2;
			$tar4 = join ",", sort keys %tar4;
		}else{
			$tad = ".", $tar1 = "."; $tar2 = "."; $tar3 = "."; $tar4 = ".";
		}
		$tar1 = "." if !$tar1; $tar2 = "." if !$tar2; $tar3 = "." if !$tar3; $tar4 = "." if !$tar4;
		print O1 "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$exon\t$noncoding\t$tar1\t$tad\n";
		print O2 "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$exon\t$noncoding\t$tar2\t$tad\n";
		print O3 "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$exon\t$noncoding\t$tar3\t$tad\n";
		print O4 "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$exon\t$noncoding\t$tar4\t$tad\n";
	}
	close I || die;
	close O1 || die;
	close O2 || die;
	close O3 || die;
	close O4 || die;
	return ();
}




1;
