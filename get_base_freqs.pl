#!/usr/bin/perl -w

use Bio::SeqIO;

my $file = shift;
my $seqin = Bio::SeqIO->new(-file => $file);
while ( my $seqobj = $seqin->next_seq() ) {
  my $seq = $seqobj->seq();
  my $id = $seqobj->id;
  #print "$id\t$seq\n";
  my $a = ($seq =~ tr/[Aa]//);
  my $c = ($seq =~ tr/[Cc]//);
  my $g = ($seq =~ tr/[Gg]//);
  my $t = ($seq =~ tr/[Tt]//);

  my $total = $a+$c+$g+$t;

  my $a_freq = $a/$total;
  my $c_freq = $c/$total;
  my $g_freq = $g/$total;
  my $t_freq = $t/$total;

  print join("\t", $id, $a_freq, $c_freq, $g_freq, $t_freq), "\n";
}
