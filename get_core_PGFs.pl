#! /usr/bin/env perl 
use strict;
use Data::Dumper;
use Getopt::Long;
use gjostat;

my $usage = 'get_core_PGFs.pl <fams_file
             /vol/patric3/fams/2016-0904/renumbered2.1.1 is an example of a fams file.
            
             -d output dir, default = "Genus_Core_PGFs";
             -b optional file of bad genome identifiers to avoid (newline delimited)
			     
			     creates a directory called "Genus_Core_PGFs" which contains a file
			     named for each genus. Each tab-delimited file contains:
			      PGF Genus AVG_PROT_Len STD_PROT_LEN ANNO
             ';
             
my ($help,  $badfile);
my $out_dir = "Genus_Core_PGFs";
my $opts = GetOptions('b=s' => \$badfile,
					  'd=s' => \$out_dir,
                      'h'   => \$help);

if ($help){die "$usage\n";}

mkdir $out_dir;

# read in a bad file if it exists
my %bad;
if ($badfile)
{
	open (IN, "<$badfile"), or die "could not open bad file\n";
	%bad = map{chomp; $_, 0}(<IN>);
	close IN;
}


my $previous = 0;
my %gene_count;
my $pgf_genome_count = {};
my $fam_data = {};

while (<>)
{
	chomp;
	my ($pgf, $n1, $n2, $peg, $len, $ann, $plf, $genus, $plf2) = split /\t/;
	
	my $gid = $peg;
	$gid =~ s/^fig\|//g;
	$gid =~ s/\.peg.+//g;
	
	$pgf =~ s/GF/PGF_/g;

	unless(exists $bad{$gid})
	{
		
		if ($genus eq $previous)
		{
			$pgf_genome_count->{$pgf}->{$gid} ++;
			$gene_count{$gid} ++;
			push @{$fam_data->{$pgf}->{LENGTHS}}, $len; 
			$fam_data->{$pgf}->{ANNO} = $ann; 
		}

		elsif ($genus ne $previous)
		{
			if ($previous)
			{
				# print the results (if they exist)
				process_fams($previous, $pgf_genome_count, $fam_data, \%gene_count);
				
				# zero out the hashes	
				$pgf_genome_count = {};
				%gene_count = ();
				$fam_data = {};
			}
			# add the current line to the hashes
			$pgf_genome_count->{$pgf}->{$gid} ++;
			$gene_count{$gid} ++;
		
			push @{$fam_data->{$pgf}->{LENGTHS}}, $len; 
			$fam_data->{$pgf}->{ANNO} = $ann; 
					
			# update $previous
			$previous = $genus;
		}
	}
}

#do the last process fams
process_fams($previous, $pgf_genome_count, $fam_data, \%gene_count);


#---------------------------------------------------------------------
#
#	process_fams ($previous, $pgf_genome_count, \%gene_count
#
#---------------------------------------------------------------------

sub process_fams
{
	my ($genus, $pgf_genome_count, $fam_data, $gene_countR) = @_;
	
	print STDERR "\n\n##PROCESSING $genus\n";
	
	# Get allowable min and max genes per genome
	my %gene_count = %$gene_countR;
	my @ngenes = values %gene_count;
	my ($avg, $std) = &gjostat::mean_stddev(@ngenes);
	my $max = ($avg + ($std * 2));
	my $min = ($avg - ($std * 2));

	# Get the usable genome list
	my %genome_list;
	foreach (keys %gene_count)
	{
		if (($gene_count{$_} > $min) && ($gene_count{$_} < $max))
		{
			$genome_list{$_} = 1;
		}
		else
		{
			print STDERR "Excluded $_\t$genus $_\tN genes = $gene_count{$_}\tAVG = $avg\tSD=$std\n";
		}
	}
	
	my $ngenomes = keys %genome_list;
	# if there are more than 4 genomes in the genus
	if ($ngenomes >= 4)
	{
		print STDERR "Processing $genus\t$ngenomes genomes\n";
		my %usable_fams;
		my $frac = 0.98;
		
		while (((keys %usable_fams) < (0.1 * $avg)) && ($frac >= 0.75))
		{
			print STDERR "$genus \tFRAC= $frac\tAVG_genes = $avg\n";
			foreach (keys %$pgf_genome_count)
			{
				my $fam = $_;				
				my $genomes_in_fam = 0;
				my $npegs = 0;
				foreach (keys %genome_list)  # the list of usable genomes
				{
					if (exists $pgf_genome_count->{$fam}->{$_})
					{
						$genomes_in_fam ++;
						$npegs += $pgf_genome_count->{$fam}->{$_};
					}
				}	
				my $frac_reps = ($genomes_in_fam/$ngenomes);

				my $frac_dupes = 0;
				unless ($npegs == 0)
				{
					$frac_dupes = ($npegs/$genomes_in_fam);
				}
				
				if (($frac_reps >= $frac) && ($frac_dupes <= 1.05))
				{
					$usable_fams{$fam} = 1; 
				}
			}
			$frac -= 0.01;
		}

		if (%usable_fams)
		{
			open (OUT, ">$out_dir/$genus");	
			foreach (sort keys %usable_fams)
			{
				my $lengthsR = $fam_data->{$_}->{LENGTHS};
				my @lengths = @$lengthsR;
				my ($avg_len, $std_len) = &gjostat::mean_stddev(@lengths);
				my $anno = $fam_data->{$_}->{ANNO};
				print  OUT  "$_\t$genus\t";
				printf OUT ("%.2f\t", $avg_len);
				printf OUT ("%.2f\t", $std_len); 
				print  OUT "$anno\n";
			}
		}
	}
}


			






























































