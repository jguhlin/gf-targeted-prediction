use strict;

my $input = $ARGV[0];
my $a;

open (my $fh, "<$input");

my $group_id = 0;
while (<$fh>) {
	chomp;
	my $outname = "out." . $group_id;
	my $output = $outname . ".fasta";
	my $out_aln = $outname . ".aln";
	my $out_prepared = $outname . ".prepared";
	my $out_prepared_aln = $outname . ".prepared.aln";
	my $out_prfl = $outname . ".prfl";

	system("rm $output");
	my @genes = split /\s+/;
	foreach my $gene (@genes) {
		$a = `samtools faidx spermatogenesis_sequences.fa "$gene" >> $output `;
	}
	$a = `mafft --quiet --maxiterate 1000 --localpair --reorder $output > $out_aln`;
	$a = `prepareAlign < $out_aln > $out_prepared 2> pa.err`;
	#$a = `mafft --quiet --maxiterate 1000 --localpair --reorder $out_prepared > $out_prepared_aln`;
	$group_id++;
	$a = `/Volumes/archive/deardenlab/guhlin/software/augustus-3.3.1/scripts/msa2prfl.pl $out_prepared > $out_prfl`;
	if (-s $out_prfl) {
		print $out_prfl . "\n";
	} else {
		$a = `rm $out_prfl`;
	}
}
