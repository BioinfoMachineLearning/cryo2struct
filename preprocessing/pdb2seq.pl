#!/usr/bin/perl -w
#################################################################
#Extract sequence from ATOM records of pdb file
#Input: pdb file
#Output: sequence in the pdb file
#Author: Jianlin Cheng
#Date: 12/30/2005
#################################################################

##############standard Amino Acids (3 letter <-> 1 letter)#######
%amino=();
$amino{"ALA"} = 'A';
$amino{"CYS"} = 'C';
$amino{"ASP"} = 'D';
$amino{"GLU"} = 'E';
$amino{"PHE"} = 'F';
$amino{"GLY"} = 'G';
$amino{"HIS"} = 'H';
$amino{"ILE"} = 'I';
$amino{"LYS"} = 'K';
$amino{"LEU"} = 'L';
$amino{"MET"} = 'M';
$amino{"ASN"} = 'N';
$amino{"PRO"} = 'P';
$amino{"GLN"} = 'Q';
$amino{"ARG"} = 'R';
$amino{"SER"} = 'S';
$amino{"THR"} = 'T';
$amino{"VAL"} = 'V';
$amino{"TRP"} = 'W';
$amino{"TYR"} = 'Y';
###################################################################

#read sequence from atom file
sub get_seq_from_atom
{
	#assume the atom file exists
	my $file = $_[0];
	open(ATOM, $file) || die "can't read atom file: $file\n";
	my @atoms = <ATOM>;
	close ATOM; 
	my $prev = -1;
	my $seq = ""; 
	while (@atoms)
	{
		my $text = shift @atoms;
		if ($text =~ /^ATOM/)
		{
			#get aa name
			#get position
			my $res = substr($text, 17, 3); 
			$res = uc($res); 
			$res =~ s/\s+//g;
			my $pos = substr($text, 22, 4);

			#if 3-letter, convert it to 1-letter.
			if (length($res) == 3)
			{
				if (exists($amino{$res}) )
				{
					$res = $amino{$res}; 
				}
				else
				{
					$res = "X"; 
					print "$file: residue is unknown, shouldn't happen.\n";
				}
			}
			if ($pos != $prev)
			{
				$seq .= $res; 
				$prev = $pos; 
			}
		}
	}
	return $seq; 
}

if (@ARGV != 1)
{
	die "need one parameters: pdb file.\n"; 
}

$pdb_file = shift @ARGV;

$seq = &get_seq_from_atom($pdb_file);
print "$seq\n";

