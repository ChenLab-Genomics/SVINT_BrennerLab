
##! /usr/bin/perl
# the above line should be adjusted according to your actual path for perl, e.g. /usr/bin/perl
use warnings;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require SV_SUB;

print "*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*\n";

$sv_file = $ARGV[0];
$file_type = $ARGV[1]; die "INPUT should be: SVfile and type-of-input (\"TXT\" or \"VCF\")\n" if $file_type !~ /^TXT$|^VCF$/;
$fi_cut = $ARGV[2];
if(!$fi_cut){
$fi_cut = 0.01;
print "Missing customized FI cut-off: use $fi_cut as default.\n";
}
print "Outputs: VCF files without FI cut-offs will be generated. TXT files will be generated at FI <= $fi_cut.\n";

$tad_bd_file = "data/tad/boundary_regions_hg19.bed"; #sorted
$tad_bd_gn_file = "data/tad/cell_bd_nearby_genes.txt"; # for each boundary in each cell, the genes between this boundary and the nearest boundaries in the same cell (or start or end point of a chromosome)

my $gene_anno;
opendir my $dir, "data/gene_anno" or die;
my @annos = readdir $dir;
for $an(@annos){
$gene_anno = "data/gene_anno/".$an if $an =~ "Homo_sapiens.GRCh37.*.PC.genes.bed\$"; # sorted
$exon_anno = "data/gene_anno/".$an if $an =~ "Homo_sapiens.GRCh37.*.PC.exons.bed\$"; # sorted
$gn2ex_anno = "data/gene_anno/".$an if $an =~ "Homo_sapiens.GRCh37.*.genes2exons.txt";
}
$re_anno = "data/reg_anno/combined_regulatory_elements_for_SVint.bed"; #sorted
$rg_all = "data/target/consolidated_RE2target_PCgenes.txt";
$kg_fi_file = "data/kg/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.hg19.vcf_based_FI.txt";
$gnomad_fi_file = "data/gnomad/gnomad_v2_sv.sites.vcf_based_FI.txt";

print "##########################################################\nSVINT v1.4 starts...\n";
print "please note: SVINT is using $kg_fi_file and $gnomad_fi_file\n";

print "Read in data...\n";
undef $same_a; undef $diff_a; undef $lines_a;
if($file_type eq "VCF"){
$sv_file_new = $sv_file;
($same_a, $diff_a, $lines_a) = SV_SUB::read_sv($sv_file);
}else{
SV_SUB::write_txt2vcf($sv_file);
$sv_file_new = $sv_file.".fake.vcf";
($same_a, $diff_a, $lines_a) = SV_SUB::read_sv($sv_file_new);
}

# same chr, diff chr; in this tool only SVs that don't cross different chr are used for analysis
undef %same_a; undef %diff_a; undef %lines_a;
%same_a = %{$same_a}; %diff_a = %{$diff_a}; %lines_a = %{$lines_a};
# $same_a{$var}{$type}=$gt; $diff_a{"$loc1\t$loc2"}{$type}=$gt; $lines_a{$id}{$line}=$var or "$loc1\t$loc2";

# filter: use all SVs
print "Annotate by gene annotation...\n";
&SV_SUB::write_bed($same_a, "$sv_file\_temp"); # output $sv_file\_temp.sorted.bed
undef $overlap_gene; undef %overlap_gene; undef %overlap_exon; undef $overlap_exon;
$overlap_gene = SV_SUB::overlap2("$sv_file\_temp.sorted.bed", $gene_anno);
%overlap_gene = %{$overlap_gene};
$overlap_exon = SV_SUB::overlap2("$sv_file\_temp.sorted.bed", $exon_anno);
%overlap_exon = %{$overlap_exon};
my %see; map{$see{$_}=1}(keys %overlap_exon);
my %to_use;
for $var(keys %same_a){
	for $tp(keys %{$same_a{$var}}){
		# next if exists $see{$var}; 
		$to_use{$var}{$tp}=$same_a{$var}{$tp};
	}
}
&SV_SUB::write_bed(\%to_use, "$sv_file\_temp2"); # output $sv_file\_temp2.sorted.bed

# annotate: with TAD boundary regions
print "Calculate SV impact on TAD boundary regions ...\n";
undef $sv2tad_bd; 
($sv2tad_bd) = SV_SUB::anno_tad1("$sv_file\_temp2.sorted.bed", $tad_bd_file, $tad_bd_gn_file);

# annotate: with regulatory elements
print "Annotate with regulatory elements ...\n";
undef $sv2re;
$sv2re = SV_SUB::overlap("$sv_file\_temp2.sorted.bed", $re_anno); #$sv2re{$var}{$re} = overlapped length ratio in RE; $var & $re: "$chr\t$stt\t$end"

# annotate with each regulatory element - target gene file:
print "Annotate with RE-target information, and their TAD boundary impact ...\n";
undef $sv2target; undef $loc_gn;
($sv2target, $loc_gn) = SV_SUB::rg_anno($sv2re, $re_anno, $rg_all, \%overlap_gene); # distal target genes for intergenic SVs, while located gene(s) as target genes for genic noncoding SVs
&SV_SUB::write_bed3($loc_gn, $gene_anno, "$sv_file\_temp3"); # output $sv_file\_temp3.sorted.bed
undef $lg2tad_bd; undef %lg2tad_bd;
($lg2tad_bd) = SV_SUB::anno_tad2("$sv_file\_RegulatoryRegion--Gene_crossTAD-cells.txt", "$sv_file\_temp3.sorted.bed", $tad_bd_file);
%lg2tad_bd = %{$lg2tad_bd};

# annotate with FI
print "Annotate with Frequency of Impact for REs calculated through 1000Genome Project SV data ...\n";
undef $fi1; undef $fi2; undef $sv2fi_tad;
$fi1 = SV_SUB::read_fi($kg_fi_file, $sv2re);
$fi2 = SV_SUB::read_fi($gnomad_fi_file, $sv2re);
$sv2fi_tad = SV_SUB::anno_fi_tad($sv2target, \%to_use, $fi1, $fi2, $lg2tad_bd);

# write out new VCF
print "Write out the new VCF and TXT files...\n";
open I, $gn2ex_anno; 
undef %gn2ex;
while(<I>){
chomp; ($gn, $ex) = split /\t/;
$gn2ex{$gn}{$ex}=1;
}
close I;

&SV_SUB::write_svint_anno_incluEx($sv_file."_SVINT_v1.4_anno_incluEx.vcf", $sv_file_new, $lines_a, $same_a, \%to_use, \%overlap_gene, \%overlap_exon, \%gn2ex, $sv2tad_bd, $sv2fi_tad);
&SV_SUB::write_txt_anno($sv_file."_SVINT_v1.4_anno_incluEx.vcf", $fi_cut, $sv_file."_SVINT_v1.4_anno_incluEx_FI$fi_cut\_impact100_2methods.txt", $sv_file."_SVINT_v1.4_anno_incluEx_FI$fi_cut\_impact100_1method.txt",$sv_file."_SVINT_v1.4_anno_incluEx_FI$fi_cut\_impact50_2methods.txt", $sv_file."_SVINT_v1.4_anno_incluEx_FI$fi_cut\_impact50_1method.txt");
print "SVINT done.\n################\n";

system "rm $sv_file\_temp.*bed"; system "rm $sv_file\_temp2.*bed"; system "rm $sv_file\_temp3.*bed*";

print "*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*\n";

