#! /usr/bin/env perl
use strict;
use gjoseqlib;
use gjostat;
use Data::Dumper;
use Getopt::Long;

my $usage = 'pattyfam_genome_quality.pl [options] >qualty_data;
			
			
			Options
			-c Genus_Core_PGFs directory [required]
			-d pattyfam data dir for place_proteins_into_pattyfams e.g., (/vol/patric3/fams/2017-0701/kmers)
			-u pattyfam url for place_proteins_into_pattyfams  d = http://spruce:6100
			-g genus [required]
			-p protein fasta file [required]
			-n no header
			-o file name for enumerated output [d = Genome.quality_data] 
			
			
			This program reads a fasta file of proteins and the file of core global Pattyfams (PGFs)
			and returns quality statistics based on the set of core PGFs that are expected. 
			The directory of core PGFs is built using get_core_PGFs.pl.
			
			
			Columns of standard output are:
			1.  Core_PGFs (total core PGFs)
			2.  Core_PGFs_Found  (core PGFs found in your genome)
			3.  Frac_Core_PGFs (fraction of core PGFs found in your genome)
			4.  Frac_Duplicated_PGFs (fraction of core PGFs that are duplicated in your genome)
				    **note that  single occurrance is not a prequesite of being a core PGF, 
				    so some duplicates may naturally occurr in the core set and thus exist your genome.
			5.  Frac_Missing_PGFs (Fraction of core PGFs that were missing)
			6.  Frac_Too_Short_PGFs (Fraction of core PGFs that are too short P-value < 0.01)
			7.  Frac_Too_Long_PGFs  (Fraction of core PGFs that are too long P-value < 0.01)
			8.  Ratio_All_PGFs/Prots (This is the fraction of proteins in your genome 
			                          that got called by PGFs, a low fraction means that you
			                          have a lot of novel or weird proteins).
		 

			The program also outputs a file with the suffix ".quality_data", which
			can be specified by the -o option. It specifically enumerates:
			
			1.  Core PGFs that are duplicated in your genome
				    Fields are: (PEG_ID, PGF_ID, Protein len, Anno)
			2.  Core PGFs that are missing in your genome
			        Fields are: (PGF_ID, Avg Len, STD Len, Anno)
			3.  Core PGFs that are too short in your genome
					Fields are: (PEG_ID, PGF_ID, Prot_Len, Avg_Len, STD_Len, Z-score, P-val, Anno)
			4.  Core PGFs that are too long in your genome
					Fields are: (PEG_ID, PGF_ID, Prot_Len, Avg_Len, STD_Len, Z-score, P-val, Anno)
			
			*****
			Note that it is critically important that the kmer build declared by the -d 
			option matches the file from which the the core PGFs were built, Otherwise PGF
			numbers will not match.
			



';

my $url = "http://spruce:6100"; 		
my ($genus, $core_dir, $fasta_file, $dataD, $help, $quality_file, $no_header);
my $opts = GetOptions('d=s'   => \$dataD,
                      'g=s'   => \$genus,
                      'c=s'   => \$core_dir,
                      'p=s'   => \$fasta_file,
                      'o=s'   => \$quality_file,
                      'h'     => \$help,
                      'n'     => \$no_header,
                      'u=s'   => \$url);

if ($help){die "$usage\n";}
unless ($genus){die "must declear a -g genus\n";}
unless ($core_dir){die "must declare a directory of core PGFs with -c\n";}
unless ($fasta_file){die "must declare a fasta file of proteins with -p\n";}
unless ($dataD){die "must declare a data directory for place_proteins_into_pattyfams with -d, e.g., /vol/patric3/fams/2017-0701/kmers \n";}


unless ($quality_file)
{
	my $genome = $fasta_file;
	$genome =~ s/(.+)(\/)(.+)$/$3/g; 
	$quality_file = "$genome.quality_data";
}


#load genus core pgfs
my $gc_file = $genus;
$gc_file =~ s/ /\\ /g;
my $pgf_stats = {};
my %core_pgfs;
open (IN, "<$core_dir/$gc_file") or die "cannot open genus core pattyfams file $gc_file\n";
while (<IN>)
{
	chomp;
	my ($pgf, $genus, $avg_len, $std_len, $anno) = split /\t/;
	$core_pgfs{$pgf} = 0;
	$pgf_stats->{$pgf}->{AVG} = $avg_len;
	$pgf_stats->{$pgf}->{STD} = $std_len;
	$pgf_stats->{$pgf}->{ANNO} = $anno;
}
close IN;
my $total_core = keys %$pgf_stats;


#get gene lengths from the fasta file
open (IN, "<$fasta_file") or die "could not open protein fasta file from -p\n";
my @seqs = &gjoseqlib::read_fasta(\*IN);
close IN;
my %prot_len;
for my $i (0..$#seqs)
{
	$prot_len{$seqs[$i][0]} = length $seqs[$i][2];
}
my $nprots = keys %prot_len;


#open the pattyfams to a file handle.
my @too_short;
my @too_long;
my @duplicated;
my %core_ids;
my %core_found;

my $total_PGFs = 0;
my $total_PLFs = 0;
open (IN, "place_proteins_into_pattyfams \-\-genus $genus $dataD $url $fasta_file |"), or die "failed to run place_proteins_into_pattyfams\n";
while (<IN>)
{
	chomp;
	if ($_ =~ /PGF/)
	{
		$total_PGFs ++;
		my ($id, $pgf, $n, $anno) = split /\t/; 
		my @pgf_data = ($id, $pgf, $prot_len{$id}, $anno);
		
		if (exists $core_pgfs{$pgf})
		{
			
			push @{$core_found{$pgf}}, \@pgf_data;
			
			my $id_anno = "$id\t$anno";
			push @{$core_ids{$pgf}}, $id_anno;
						
			#determine if it is too long or too short:
			my $avg = $pgf_stats->{$pgf}->{AVG};
			my $std = $pgf_stats->{$pgf}->{STD};
			my $len = $prot_len{$id};
			my $z = (($len - $avg)/($std + 0.000000001));
			
			if ($z < 0)
			{
				my $p = &gjostat::std_normal_le_z($z);
				if ($p < 0.01)
				{
					my @array = ($id, $pgf, $len, $avg, $std, $z, $p, $anno);
					push @too_short, \@array;
				}
			}
			if ($z > 0)
			{
				my $p = &gjostat::std_normal_ge_z($z);
				if ($p < 0.01)
				{
					my @array = ($id, $pgf, $len, $avg, $std, $z, $p, $anno);
					push @too_long, \@array;
				}
			}
		}		
	}
}
my $n_core_found = keys %core_found;

# round up duplicated pgfs
my @dups;
foreach (keys %core_found)
{
	my @array = @{$core_found{$_}};
	
	my $size = scalar @array;
	if ($size > 1)
	{
		push @dups, \@array;
	}
}
my $n_dups = scalar @dups;


# round up missing pgfs
my @missing;
foreach (keys %core_pgfs)
{
	my $pgf = $_;
	unless (exists $core_found{$pgf})
	{
		my @array = ($pgf, $pgf_stats->{$pgf}->{AVG}, $pgf_stats->{$pgf}->{STD}, $pgf_stats->{$pgf}->{ANNO}); 
		push @missing, \@array;
	}
}	
my $n_missing = scalar @missing;



unless ($total_core){die "No core genes found for $genus\n";}
my $frac_core   = ($n_core_found/$total_core);
my $frac_dup    = ($n_dups/$total_core);
my $frac_miss   = ($n_missing/$total_core);
my $n_too_short = scalar @too_short;
my $frac_short  = ($n_too_short/$total_core);
my $n_too_long  = scalar @too_long;
my $frac_long   = ($n_too_long/$total_core);
my $ratio_pgf_prot = ($total_PGFs/$nprots);

unless ($no_header)
{
	print "Core_PGFs\tCore_PGFs_Found\tFrac_Core_PGFs\tFrac_Duplicated_PGFs\tFrac_Missing_PGFs\tFrac_Too_Short_PGFs\tFrac_Too_Long_PGFs\tRatio_All_PGFs/Prots\n"; 
}

print "$total_core\t$n_core_found\t";

printf ("%.3f", $frac_core);
print "\t";
printf ("%.3f", $frac_dup);
print "\t";
printf ("%.3f", $frac_miss);
print "\t";
printf ("%.3f", $frac_short);
print "\t";
printf ("%.3f", $frac_long);
print "\t";
printf ("%.3f", $ratio_pgf_prot);
print "\n";


open (OUT, ">$quality_file");

print OUT "## Core PGFs that are duplicated:\n";
for my $i (0..$#dups)
{
	for my $j (0..$#{$dups[$i]})
	{
		print OUT join "\t", @{$dups[$i][$j]}, "\n";
	}
}

print OUT "\n\n## Core PGFs that are missing:\n";
for my $i (0..$#missing)
{
	print OUT join "\t", @{$missing[$i]}, "\n";
}

print OUT "\n\n## Core PGFs that are too short:\n";
for my $i (0..$#too_short)
{
	print OUT join "\t", @{$too_short[$i]}, "\n";
}

print OUT "\n\n## Core PGFs that are too long:\n";
for my $i (0..$#too_long)
{
	print OUT join "\t", @{$too_long[$i]}, "\n";
}












