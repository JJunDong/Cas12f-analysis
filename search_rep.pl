#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use constant LEN => 15; # repeat length
use constant WIN => 10; # sliding wind size
use constant ENDFL => 100; # 5' to end of aa seq containing repeat

use constant COV => 0.85; # coverage of short repeat 
use constant LOH => 3; # length of overhang
use constant NUM_N => 2; # number of N allowed
 
my $threads = 2;
my $infile  = shift;
my $outfile = shift;
my $outSeqFile = $outfile;
$outSeqFile=~s/\.txt$/\.fasta/;
my $in = Bio::SeqIO->new(-file =>$infile, -format => 'fasta');
my $fout = Bio::SeqIO->new(-file => ">$outSeqFile", -format => 'fasta');
open(REOUT, ">$outfile") or die "Cannot write $outfile: $!";
my %seqLibOut=();

print STDERR "Starting searching IR and Repeat\n";
while(my $seqObj = $in->next_seq){
	my $id = $seqObj->id;
	print STDERR "\nFor $id\n";
	my $db_file = $id.'_bln.fasta';
	my $out = Bio::SeqIO->new(-file => ">$db_file", -format => 'fasta');
	$out->write_seq($seqObj);
	
	print STDERR "\tStart Blast\n";
	my $makeblastdb_cmd = "makeblastdb -in $db_file -dbtype nucl -out tmp_db -title tmp_db";
	system($makeblastdb_cmd) == 0 or die "Fail to compile blast database";
	
	my $length = $seqObj->length;
	my $desc = $seqObj->desc;
	my @pos=$desc=~/(\d+)/g;
	# print $id,"\t", $length,"\t",join("\t", @pos),"\n";
	my $start = 1;
	my $end = $start + LEN - 1;
	
	my $bln_output = 'blastnSearchOut_'.$id.'.txt';
	open(OUT, " | blastn -db tmp_db -num_threads $threads -task 'blastn-short' -outfmt '6 qseqid qstart qend sseqid sstart send length nident evalue bitscore qframe sframe' > $bln_output");
	while($end < $length){
		unless($start > $pos[0] and $end < $pos[1] - ENDFL){
			my $subseq = $seqObj->subseq($start, $end);
			my $tmp_id = $id.'_subseq:'.$start.'-'.$end;
			my $num_n = $subseq =~tr/Nn/NN/;
			print OUT ">$tmp_id\n",$subseq,"\n" if($num_n <= NUM_N);
		}
		$start+=WIN - 1;
		$end = $start + LEN - 1;
	}
	close(OUT);	
	unlink ('tmp_db.nhr', 'tmp_db.nin', 'tmp_db.nsq') or die "failed to remove tmp_db files"; 
	
	## parse bln output
	print STDERR "\tParsing Blast\n";
	open(IN, $bln_output) or die "Cannot open $bln_output: $!";
	my @content = <IN>;
	close(IN);
	chomp @content;
	foreach my $bln_str (@content){
		my @bln_str = split/\s+/, $bln_str;
		my ($qid, $qstart, $qend, $sid, $sstart, $send, $aln_length, $nident, $evalue, $score, $qframe, $sframe) = @bln_str;
		my ($q_start, $q_end) = $qid=~/subseq:(\d+)-(\d+)$/;
		# print "$q_start, $q_end\n";
		next if($q_start == $sstart and $q_end == $send); # remove self-bln record
		if($nident / LEN >= COV and LEN - ($qend - $qstart +1) <= LOH){
			# panlindrome
			if($q_end < $pos[0] and $sstart >= $pos[1] - ENDFL and $qframe ne $sframe){
				# $qid = (split/_/,$qid)[0];
				print REOUT "$qid\t$q_start\t$q_end\t", $seqObj->subseq($q_start, $q_end),"\t", 
				            join("\t", @bln_str[3..5]),"\t", $seqObj->subseq($send,$sstart),"\tIR\n";
				$seqLibOut{$id}=$seqObj;
			}elsif($q_start >= $pos[1] - ENDFL and $sstart > $q_end and $qframe eq $sframe){
				# $qid = (split/_/,$qid)[0];
				print REOUT "$qid\t$q_start\t$q_end\t", $seqObj->subseq($q_start, $q_end),"\t",
				            join("\t", @bln_str[3..5]),"\t", $seqObj->subseq($sstart,  $send),"\tRepeat\n";
				$seqLibOut{$id}=$seqObj;
			}
		}
	}
	print STDERR "\tDone for $id\n";
	
}
close(REOUT);

print STDERR "\n Writing results\n";
foreach my $id (sort keys %seqLibOut){
	my $seqObj = $seqLibOut{$id};
	$fout->write_seq($seqObj);
}

# opendir(PATH, '.') or die "Cannot open \.:$!";
# my @bln_files = grep(/^blastnSearchOut/, readdir PATH);
# foreach my $bln_file (@bln_files){
# 	open(IN, $bln_file) or die "Cannot open $bln_file: $!";
# 	my @content = <IN>;
# 	close(IN);
# 	chomp @content;
# 	foreach (@content){
#
# 	}
#
# }





















