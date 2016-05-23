#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Getopt::Long;
use Pod::Usage; 
 
my $man;
my $help;
my $gff;
my $prefix;
my $fasta; 
# /homedir/droc/GenomeHub/cc2-login/converter/gff2fasta.pl --gff /homedir/droc/GenomeHub/sample/Sh_205G05.gff3 --fasta /homedir/droc/GenomeHub/sample/Sh_205G05.fna --prefix Sh_205G05
my $usage = q/
 gff2fasta.pl --gff <gff_file> --fasta <fasta_file> --prefix <output_prefix>

Parameters  
	--gff               GFF3 file  [required] 
	--fasta             FASTA file [required]
	--prefix            Prefix for  
	--help
/;
GetOptions(
    'gff=s'    => \$gff,
    'fasta=s'  => \$fasta,
    'prefix=s' => \$prefix,   
    'help|h|?' => \$help
)  or pod2usage(-message => $usage);
if ($help) { pod2usage(-message => $usage); } 
 
if ($gff eq "") {
    print "Warn :: --gff is empty\nPlease specify the GFF file\n";
    print $usage;
    exit 0;
}
if ($fasta eq "") {
    print "Warn :: --fasta is empty\nPlease specify the FASTA file\n";
    print $usage;
    exit 0;
}
my $output_gene = $prefix ."_gene.fna";
my $output_cds = $prefix ."_cds.fna";
my $output_cdna = $prefix ."_cdna.fna";
my $output_protein = $prefix ."_protein.faa";

# Input file
my $fasta_in = new Bio::SeqIO(
	-file => $fasta,
	-format => 'fasta'
);
my $gff_in = new Bio::Tools::GFF(
	-file => $gff,
	-gff_version => 3
);
my %is_in_list;
if ($list_gene) {
	open(IN,$list_gene);
	while(<IN>){
		chomp;
		$is_in_list{$_} = 1; 
	}
	close IN;
}
if ($list_cds){
	open(IN,$list_cds);
	while(<IN>){
		chomp;
		$is_in_list{$_} = 1; 
	}
	close IN;
}  
# Output file
###### Output type description ######
# gene 			- the entire gene sequence (including UTRs and introns)
# cds 			- translated sequence (starting with ATG and ending with a stop codon included)
# cdna 			- transcribed sequence (devoid of introns, but containing untranslated exons)
# protein 		- cds translated (includes a * as the stop codon)

my $gene_out = new Bio::SeqIO(
	-file => ">$output_gene",
	-format => 'fasta'
);
my $cds_out = new Bio::SeqIO(
	-file => ">$output_cds",
	-format => 'fasta'
);
my $cdna_out = new Bio::SeqIO(
	-file => ">$output_cdna",
	-format => 'fasta'
);
my $protein_out = new Bio::SeqIO(
	-file => ">$output_protein",
	-format => 'fasta'
);

# Read FASTA
 
my %fasta;
while (my $seqobj = $fasta_in->next_seq) {
	$fasta{$seqobj->display_id} = $seqobj;
}
$fasta_in->close;


# Read & parse GFF3
my %gene;
my %exon;
my %cds;
my %desc;
while(my $feature = $gff_in->next_feature) {
	if ($feature->primary_tag() eq "gene") {
        my ($gene_id) = $feature->get_tag_values("ID");
		if ($list_gene) { 
            if (defined $is_in_list{$gene_id}) {
				push @{$gene{$feature->seq_id}} , $feature;
			} 
        }
        else {
			push @{$gene{$feature->seq_id}} , $feature;
		}
    }
	elsif ($feature->primary_tag() eq "mRNA") {
        my ($mrna_id) = $feature->get_tag_values("ID");
        ($note) = $feature->get_tag_values("Note") if $feature->has_tag("Note");
		$desc{$mrna_id} = $note if $note;
	}
	elsif ($feature->primary_tag() eq "CDS") {
		my ($mrna_id) = $feature->get_tag_values("Parent");
		if ($list_cds) {
            if (defined $is_in_list{$mrna_id}) {
				push @{$cds{$feature->seq_id}{$mrna_id}} , $feature;
			}
        }
        else {
			push @{$cds{$feature->seq_id}{$mrna_id}} , $feature;
		}
	} 
	elsif ($feature->primary_tag() eq "exon") {
		my ($mrna_id) = $feature->get_tag_values("Parent");
		if ($list_cds) {
            if (defined $is_in_list{$mrna_id}) {
				push @{$exon{$feature->seq_id}{$mrna_id}} , $feature;
			}
        }
        else {
			push @{$exon{$feature->seq_id}{$mrna_id}} , $feature;
		}
	}
}
$gff_in->close;

# Write Fasta
my $cpt_gene    = 0;
my $cpt_cds     = 0;
my $start_codon = 0;
my $stop_codon  = 0;
my @bad_sequence;
foreach my $seq_id (keys %fasta) {
	my $seqobj = $fasta{$seq_id};
	my $strand;
	foreach my $feature (sort {$a->start <=> $b->start} @{$gene{$seq_id}}){
		$cpt_gene++;
        my ($gene_id) = $feature->get_tag_values("ID");
        ($note) = $feature->get_tag_values("Note") if $feature->has_tag("Note");
		my $seqobj_gene = $seqobj->trunc($feature->start,$feature->end);
		my $strand = $feature->strand;
        if ($strand =~ /-/) {
            $seqobj_gene = $seqobj_gene->revcom();
        }	 
		$seqobj_gene->display_id($gene_id);
		$seqobj_gene->desc($note) if $note;
		$gene_out->write_seq($seqobj_gene);
		
	}
	foreach my $mrna_id (keys %{$cds{$seq_id}}){
		$cpt_cds++;
		my $prot_id = $mrna_id;
		$prot_id =~ s/t/p/;
		my @cds = @{$cds{$seq_id}{$mrna_id}};
		my @exon = @{$exon{$seq_id}{$mrna_id}};
		my $desc = $desc{$mrna_id} if defined $desc{$mrna_id};
		my $cds_seq;
		foreach my $feature (sort{$a->start <=> $b->start} @{$cds{$seq_id}{$mrna_id}}) {	
			$cds_seq .= $seqobj->subseq($feature->start,$feature->end); 
			$strand = $feature->strand;
        } 
		my $seqobj_cds = Bio::PrimarySeq->new(
            -seq        => $cds_seq, 
            -display_id => $mrna_id
        );
        if ($strand =~ /-/) {
            $seqobj_cds = $seqobj_cds->revcom();
        }
        if ($seqobj_cds->seq =~ /^ATG.*/){
			$start_codon++;
        }
		else {
			push @bad_sequence , $mrna_id;
		}
        if ($seqobj_cds->seq =~ /.*(TAG|TAA|TGA)$/){
			$stop_codon++;
		}
		$seqobj_cds->desc($desc) if $desc;
        my $seqobj_protein = $seqobj_cds->translate();
		$seqobj_protein->display_id($prot_id);
        $cds_out->write_seq($seqobj_cds);
        $protein_out->write_seq($seqobj_protein);	
		my $cdna_seq;
		foreach my $feature (sort{$a->start <=> $b->start} @{$exon{$seq_id}{$mrna_id}}) {	
			$cdna_seq .= $seqobj->subseq($feature->start,$feature->end); 
        }	
		my $seqobj_cdna = Bio::PrimarySeq->new(
            -seq        => $cdna_seq, 
            -display_id => $mrna_id
        );
        if ($strand eq '-') {
            $seqobj_cdna = $seqobj_cdna->revcom();
        }
        $cdna_out->write_seq($seqobj_cdna);
    }
}
$gene_out->close;
$cds_out->close;
$cdna_out->close;
$protein_out->close;
print "Number of gene(s) : ". $cpt_gene ."\n";
print "Number of transcript : " . $cpt_cds ." (with ATG ". $start_codon ."; with stop_codon ". $stop_codon .")\n";
print join("\n",@bad_sequence) ,"\n";
