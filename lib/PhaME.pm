#!/usr/bin/env perl

use strict;
use warnings;
use FileHandle;
use File::Basename;
use Cwd;

package PhaME;


sub check {

    # Rune some checks
    # Check for first time run
    # Check if the previous run finished properly
    my $dir       = shift;
    my $refdir    = shift;
    my $time      = shift;
    my $data      = shift;
    my $reference = shift;
    my $log       = shift;
    my $project   = shift;
    my $list      = $dir . '/working_list.txt';
    my $wdir      = $dir . '/results';
    my $progress  = 0;
    my $alignment = 0;
    my %refcheck;

    if ( $time == 1 ) {

        my @overwrite = glob("$wdir/RAxML_*.$project\_all");
        @overwrite = glob("$wdir/RAxML_*.$project\_cds");
        if ( scalar @overwrite > 0 ) {
            print
                "*WARNING* RAxML trees with the name $project already exist. They will be overwritten.\n";
            system("rm -f $wdir/RAxML*");
        }

        if ( -e $log && !-z $log ) {
            open( LOG, "$log" ) || die "$!";
            while (<LOG>) {
                chomp;
                if (/NUCmer on all reference genomes complete/) {
                    $progress = 1;
                }
                if (/NUCmer on all contig\/draft genomes complete/) {
                    $progress = 2;
                }
                if (/Read Mapping complete/) { $progress = 3; }
                if (/SNP alignment complete/) {
                    $progress  = 4;
                    $alignment = 1;
                }
                if ( /Tree phylogeny complete/ && $alignment == 1 ) {
                    $progress = 5;
                }

#         if (!/Tree phylogeny complete/ && -e "$wdir/RAxML_*.$project"){`rm $wdir/RAxML_*.$project`;}
            }
        }

        if ( -e $log && $data <= 2 ) {
            if ( $progress == 3 ) { return ($progress); }
            if ( $progress < 3 ) {
                print "\n\tWarning: Previous run not complete.\n\t";
                return ($progress);
            }
        }
        if ( -e $log && $data >= 3 && $data <= 5 ) {
            if ( $progress == 4 ) { return ($progress); }
            if ( $progress < 4 ) {
                print "\n\tWarning: Previous run not complete.\n\t";
                return ($progress);
            }
        }
        if ( -e $log && $data == 6 ) {
            if ( $progress == 5 ) { return ($progress); }
            if ( $progress < 5 ) {
                print "\n\tWarning: Previous run not complete.\n\t";
                return ($progress);
            }
        }
        if ( !-e $log ) { return ("0"); }
    }    # time==1

    # Check if all snps and gaps files are present in appropriate directories
    if ( $time == 2 ) {
        my $pref_name = $reference;
        $pref_name =~ s/\W/_/g;    # Replacing special characters with _

        #print "This is the 2nd time around!\n";
        my $snpdir = "$wdir/snps";
        my $gapdir = "$wdir/gaps";
        if ( !-e $list ) {
            print
                "\n$list not found.\nChange control file to run the entire pipeline.\n";
            return ("1");
        }
        else {
            open( LIST, $list );
            while (<LIST>) {
                chomp;

                my $snpsfile = "$snpdir/$pref_name\_$_.snps";
                my $gapsfile = "$gapdir/$pref_name\_$_.gaps";
                if (/contig/) {
                    $snpsfile = "$snpdir/$pref_name\_$_.snps";
                    $gapsfile = "$gapdir/$pref_name\_$_.gaps";
                }
                if (/contigs/) {
                    $snpsfile = "$snpdir/$pref_name\_$_" . "s.snps";
                    $gapsfile = "$gapdir/$pref_name\_$_" . "s.gaps";
                }
                if (/(.+)_p|sread/) {
                    $snpsfile = "$snpdir/$pref_name\_$1.vcf";
                    $gapsfile
                        = glob "$gapdir/$pref_name\_$1\_$pref_name*.gaps";
                }

                if ( $pref_name ne $_ && !-e $snpsfile ) {
                    print "$snpsfile\n";
                    print
                        "\nsnps directory not complete.\nChange control file to finish previous run.\n";
                    return ("1");
                }
                if ( $pref_name ne $_ && !-e $gapsfile ) {
                    print "$gapsfile\n";
                    print
                        "\ngaps directory not complete.\nChange control file to finish previous run.\n";
                    return ("1");
                }
            }
            return ("0");
        }
        close LIST;
    }    # time==2
}    # sub check

# Concatenate multiple reference chromosomes from one genome into one mega chromosome
# Saves them in a separate directory "/files" under the working directory
# Returns the name of the genome file and the total length of the genome
sub prepareComplete {
    my $wdir = shift;
    my $file = shift;
    my $name = shift;
    my $sequence;
    $name =~ s/\W/_/g;

    #print "$name\t$wdir/files/$name.fna\n";
    `mkdir -p $wdir/files`;
    open( OUT, ">$wdir/files/$name.fna" ) || die "$!";
    print OUT ">$name\n";
    open( IN, $file ) || die "$!";
    while (<IN>) {
        chomp;
        $_ =~ s/\r\n|\n|\r/\n/g; #this line
        if ( !/^>/ ) { $sequence .= $_; }    #print OUT $_;}
    }
    print OUT "$sequence\n";
    close OUT;
    close IN;
    return $name, length $sequence;
}    # prepareComplete

sub prepareContig {

    # Change contigs names
    # Saves them in a separate directory "\files" under the working directory
    # Returns the name of the contig file and the number of contigs present
    my $dir  = shift;
    my $file = shift;
    my $name = shift;
    my ( $header, @seq );
    my $sequence;
    my $count = 1;
    $name =~ s/\W/_/g;
    my $contig     = $name . '_contig';
    my $outfile    = $dir . '/files/' . $name . '_contig.fna';
    my $total_size = 0;

    #print "$contig\t$outfile\n";

    open( OUT, ">$outfile" ) || die "$!";
    my $fh = FileHandle->new($file) || die "$!";
    if ( $fh->open("<$file") ) {
        $/ = ">";
        while (<$fh>) {
            $_ =~ s/\>//g;
            $_ =~ s/\r\n|\n|\r/\n/g; # this line
            unless ($_) { next; }
            ( $header, @seq ) = split /\n/, $_;
            $sequence = join "", @seq;
            $header = $name . '_' . $count;

            #      print ">$header\n$sequence\n";
            print OUT ">$header\n$sequence\n";
            $count++;
            $total_size += length $sequence;
        }
        $/ = "\n";
        $fh->close;
        close OUT;
    }
    return $contig, $total_size;
}

# Identifies gap coords in reference genomes
# Gaps identified in NUCmer
sub identifyGaps {
    my $dir       = shift;
    my $list      = shift;
    my $name      = shift;
    my $type      = shift;
    my $project   = shift;
    my $error     = shift;
    my $log       = shift;
    my $gapdir    = $dir . '/gaps';
    my $repeatdir = $dir . '/stats';
    my %query;
    my $line = 0;
    my $all_gapfile;
    my $gap_start;
    my $gap_end;

    open( OUT, ">>$log" );

    if ( $type =~ /map/ ) {
        $all_gapfile = "$dir\/$project\_mapping_gaps.txt";
    }
    elsif ( $type =~ /snp/ ) {
        $all_gapfile = "$dir\/$project\_all_gaps.txt";
    }
    open( GAP, ">$all_gapfile" ) || die "$!";

    open( LIST, "$list" ) || die "$!";
    while (<LIST>) {
        chomp;
        $query{$_}++;
    }
    close LIST;

    opendir( DIR, "$gapdir" );
    while ( my $gaps = readdir(DIR) ) {

    if (   $gaps =~ /^$name\_norepeats\_(.+)\_norepeats\.gaps/
        || $gaps =~ /^$name\_(.+_contig)s?\.gaps/
        || $gaps =~ /^$name\_(.+)\.gaps/ ){
        if ( exists $query{$1} ) {
            $line = 0;
            my $gapfile = "$gapdir/$gaps";
            open( IN, $gapfile ) || die "$!";
            while (<IN>) { $line++; print GAP "$_"; }
            close IN;
            if ( $line == 0 ) {
                print OUT "Empty Gap Files: $gapfile\n";
                $line = 0;
            }
        }
    }
    
    if ( $type =~ /snp/ ) {
  #      if ($gaps=~ /^$name\_norepeats\_(.+)\_?$name?\_?[\d+\_\d+]?\.gaps$/){
        if ( $gaps =~ /^$name\_(.+)\_$name(\_\d+\_\d+)?\.gaps$/ ) {
            my $query = $1;
            my @read_types = ( "_pread", "_sread", "_read" );
            foreach my $type (@read_types) {
                my $tmp = $query . $type;
                if ( exists $query{$tmp} ) {
                    my $gap_file = "$gapdir/$gaps";
                    $line = 0;
                    open( IN, $gap_file ) || die "$!";

                        # print "Read Mapping Gaps $gaps\n";
                    while (<IN>) {
                        chomp;
                        $line++;
                        next if (/Start\s+End\s+Length.+/);
                        my ( $start, $end, $length, $ref ) = split "\t",
                        $_;
                            if ( $ref =~ /$name\_(\d+)\_(\d+)$/ ) {
                                $gap_start = $start + $1 - 1;
                                $gap_end   = $1 + $end - 1;
                            }
                            else {
                                $gap_start = $start;
                                $gap_end   = $end;
                            }
                            print GAP
                                "$name\t$gap_start\t$gap_end\t$length\t${query}$type\n";
                        }
                        close IN;
                        if ( $line == 1 ) { `rm $gap_file`; $line = 0; }
                    }
                }    #foreach @read_types
            }
        }

        #   last OUTER;
    }
    closedir DIR;

    opendir( REPEAT, "$repeatdir" );
    while ( my $repeats = readdir(REPEAT) ) {
        if ( $repeats =~ /$name\_repeat_coords\.txt/ ) {

            #      print "$repeats\n";
            #   if ($gaps=~ /repeats_coords\.txt/){
            my $repeatfile = "$repeatdir/$repeats";
            open( IN, $repeatfile ) || die "$!";

            #      print "Repeats\n";
            while (<IN>) {
                chomp;
                my ( $ref, $temp, $start, $end, $length ) = split "\t", $_;
                print GAP "$ref\t$start\t$end\t$length\tREPEATS\n";
            }
            close IN;
        }

        #   last OUTER;
    }
    closedir REPEAT;
    close OUT;
    return $all_gapfile;
}
#------------------------------------------------------------------------------#
# Identify CDS coords
sub codingRegions {
    my $dir        = shift;
    my $annotation = shift;
    my $name       = shift;
    my $start;
    my $end;
    my $gap_start = 1;
    my $gap_end;
    my $source_start = 1;
    my $source_end   = 0;
    my %CDS;
    my $line;
    my $temp;

    my $outfile = $dir . "/noncoding.txt";
    my $coding  = $dir . "/CDScoords.txt";
    my $coding_GFF = $dir . "/$name\_cds.gff";

    open( OUT, ">$outfile" );
    open( CDS, ">$coding" );
    open( CDS_GFF, ">$coding_GFF" );

    my $first       = 1;
    my $permutation = 0;

    open( IN, "$annotation" ) || die "$!";
    while (<IN>) {
        chomp;
        if (/##sequence-region/) {
            $permutation = $permutation + $source_end;
            ( $line, $temp, $source_start, $source_end ) = split " ", $_;
        }
        if ( !/^#/ ) {
            my ($chr_name,  $source, $method, $start, $stop,
                $score, $strand, $phase,  $field
            ) = split "\t", $_;

            my @group = split ";", $field if ($field);

            if ( $method and $method =~ /CDS/ ) {
                $start = $start + $permutation;
                $stop  = $stop + $permutation;

                print CDS "$chr_name\t$start\t$stop\t";
                print CDS_GFF "$name\t$source\t$method\t$start\t$stop\t$score\t$strand\t$phase\t$field";
                $CDS{$start} = $stop;
                foreach (@group) {
                    if ( /product=(.+)/ || /description=(.+)/ ) {
                        print CDS "$1\t";
                    }
                    if ( /gene=(.+)/ ) {
                        print CDS "$1\t";
                    }
                }
                print CDS_GFF "\n";
                print CDS "\n";
            }
        }
    }

    my $prev = 0;
    my $last = 0;
    foreach my $begin ( sort { $a <=> $b } keys %CDS ) {
        my $end = $CDS{$begin};
        if ($first) {
            if ( $begin == 1 ) { $gap_start = $end + 1; }
            else {
                $gap_end = $begin - 1;
                if ( $gap_start < $gap_end ) {
                    print OUT "$name\t$gap_start\t$gap_end\tnoncoding\n";
                }
            }
            $first = 0;
        }
        else {
            $gap_end = $begin - 1;
            if ( $gap_start < $gap_end ) {
                print OUT "$name\t$gap_start\t$gap_end\tnoncoding\n";
            }
        }
        if ( $begin == $prev ) {
            $gap_start = $prev - 1;
            $gap_end   = $prev - 1;
            print OUT "$name\t$gap_start\t$gap_end\tnoncoding\n";
        }

        $gap_start = $end + 1;
        $prev      = $end + 2;
        $last      = $end;
    }
    if ( $last < $source_end ) {
        $gap_start = $last + 1;
        print OUT "$name\t$gap_start\t$source_end\tnoncoding\n";
    }
}

sub clean {
    my $dir = shift;
    print "\nCleaning up...\n";
    system(
        "rm -f $dir/*.pileup $dir/*.bam $dir/*.bcf $dir/*.mgaps $dir/*.ntref $dir/*.sam $dir/*.delta $dir/*.*filter"
    );
}

sub movefiles {
    # Function to move files to directories
    my $dir = shift;
    my $trees = shift;
    my $reads = shift;
    
    print "\nFinalizing...\n";
    system(
        "mkdir $dir/alignments $dir/tables $dir/miscs; echo $dir/*.fna | xargs mv -t $dir/alignments/; echo $dir/*.txt | xargs mv -t $dir/tables; echo $dir/*.delta | xargs mv -t $dir/miscs; echo $dir/*.*filter | xargs mv -t $dir/miscs"
    );
    if ($trees > 0){
        system("mkdir $dir/trees");
        if ($trees == 1 || $trees == 4){
            system("mv $dir/*.fasttree $dir/trees/");
        }
        elsif ($trees == 2 || $trees == 4){
            system("mv $dir/RAxML*.* $dir/trees/");
        }
        elsif ($trees == 3 || $trees == 4){
            system("mv $dir/*.fna.* $dir/trees/");
        }
    }
    if ($reads == 2 || $reads == 4 || $reads == 5 || $reads == 6){
        system(
            "mkdir $dir/maps; mv $dir/*.coverage $dir/maps/; mv $dir/*plots.pdf $dir/maps/;mv $dir/*.bam* $dir/maps/");
    
    }
}

#-------------------------------------------------------------------------
#  SCRIPT RUNNERS

# Run NUCmer on reference genomes
sub completeNUCmer {
    my $reference = shift;
    my $indir     = shift;
    my $bindir    = shift;
    my $list      = shift;
    my $code      = shift;
    my $thread    = shift;
    my $error     = shift;
    my $log       = shift;
    my $buildSNPdb = shift;
    my $outdir    = $indir . '/results';

    open( OUT, ">>$log" );
    print OUT "\n";
    my $nucmer
        = "runNUCmer.pl -r $reference -q $indir -d $outdir -t $thread -l $list -c $code -b $buildSNPdb 2>$error > $log\n\n";
    print OUT $nucmer;
    if ( system($nucmer) ) { die "Error running $nucmer.\n"; }

    close OUT;
}

# Run NUCmer on contigs
# Needs a reference
sub contigNUCmer {
    my $indir     = shift;
    my $bindir    = shift;
    my $list      = shift;
    my $code      = shift;
    my $thread    = shift;
    my $reference = shift;
    my $type      = shift;
    my $error     = shift;
    my $log       = shift;
    my $outdir    = $indir . '/results';

    open( OUT, ">>$log" );
    print OUT "\n";
    my $con_nucmer
        = "runContigNUCmer.pl -r $reference -q $indir -d $outdir -l $list -t $thread -y $type 2>>$error >> $log\n\n";
    print OUT $con_nucmer;
    if ( system($con_nucmer) ) { die "Error running $con_nucmer.\n"; }
    close OUT;
}

# Removes gaps using a gap file
sub removeGaps {
    my $bindir    = shift;
    my $reference = shift;
    my $readgaps  = shift;
    my $error     = shift;
    my $log       = shift;
    
    open( OUT, ">>$log" );
    print OUT "\n";
    my $remove = "removeGaps.pl $reference $readgaps\n\n";
    print OUT $remove;
    if ( system($remove) ) { die "Error running $remove.\n"; }
    close OUT;
}

# Runs bowtie on paired-end reads

sub readsMapping {
    my $indir     = shift;
    my $bindir    = shift;
    my $list      = shift;
    my $thread    = shift;
    my $name      = shift;
    my $error     = shift;
    my $aligner	  = shift;
    my $ploidy    = shift;
    my $snp_filter    = shift;
    my $log       = shift;
    my $outdir    = $indir . "/results";
    my $reference = $outdir . '/temp/' . $name . '.fna';
    my $type;

    open( OUT, ">>$log" );
    if ( !-e $reference || -z $reference ) {
        $reference = $indir . '/files/' . $name . '.fna';
    }
    print OUT "\n";
    my $map
        = "runReadsMapping.pl -r $reference -q $indir -d $outdir -t $thread -l $list -a $aligner -p $ploidy -s $snp_filter 2>>$error >> $log\n\n";
    print OUT $map;
    if ( system($map) ) { die "Error running $map.\n"; }

    opendir( CLEAN, "$outdir" );
    while ( my $file = readdir(CLEAN) ) {
        if ( $file =~ /.+\.vcf$/ ) {
            my $vcf_file = $outdir . '/' . $file;
            `mv $vcf_file $outdir/snps`;
            # print "Moved $file to the snps directory\n";
        }
    }
    closedir CLEAN;
    return ("Read Mapping complete");
    close OUT;

}

sub buildSNPDB {
    my $outdir     = shift;
    my $bindir     = shift;
    my $reference  = shift;
    my $list       = shift;
    my $project    = shift;
    my $signal     = shift;
    my $error      = shift;
    my $log        = shift;
    my $gap_cutoff = shift;

    open( OUT, ">>$log" );
    print OUT "\n";
    my $SNPdb
        = "buildSNPDB.pl -i $outdir -r $reference -l $list -p $project -c $signal -g $gap_cutoff 2>>$error >> $log\n\n";
    print OUT $SNPdb;
    if ( system($SNPdb) ) { die "Error running $SNPdb.\n"; }

    close OUT;
    return ("SNP database complete");
}
#------------------------------------------------------------------------------------------------------------------------------------------------#
sub buildTree {
    my $bindir = shift;
    my $outdir = shift;
    my $thread = shift;
    my $tree   = shift;
    my $name   = shift;
    my $bsignal = shift; # signal showing if bootstrap should be done or not.
    my $bootstrap = shift;
    my $error  = shift;
    my $log    = shift;

    open( OUT, ">>$log" );

    if ( $tree == 1 || $tree == 4 ) {
        print OUT "Reconstructing phylogeny using FastTree\n";
        my $fasttree
            = "export OMP_NUM_THREADS=$thread; FastTree -quiet -nt -gtr < $outdir/$name\_all_alignment.fna > $outdir/$name\.fasttree 2>>$error \n\n";
            # = "export OMP_NUM_THREADS=$thread; FastTree -quiet -nt -gtr < $outdir/$name\_snp_alignment.fna > $outdir/$name\.fasttree 2>>$error \n\n";
        print OUT $fasttree;
        if ( system($fasttree) ) { die "Error running $fasttree.\n"; }

    }
    if ( $tree == 2 || $tree == 4 ) {
        print OUT "Reconstructing phylogeny using RaxML\n";
        print OUT "\n";
        my $raxml
            = "raxmlHPC-PTHREADS -p 10 -T $thread -m GTRGAMMAI -s $outdir/$name\_all_alignment.fna -w $outdir -n $name 2>>$error >> $log\n\n";
            # = "raxmlHPC-PTHREADS -p 10 -T $thread -m GTRGAMMAI -s $outdir/$name\_snp_alignment.fna -w $outdir -n $name 2>>$error >> $log\n\n";
        print OUT $raxml;
        if ( system($raxml) ) { die "Error running $raxml.\n"; }
        if ( $bsignal == 1 ) {
            open( OUT, ">>$log" );
            my $bootTrees
                = "raxmlHPC-PTHREADS -p 10 -T $thread -m GTRGAMMAI -b 10000 -t $outdir/RAxML_bestTree.$name -s $outdir/$name\_all_alignment.fna -w $outdir -N $bootstrap -n $name\_b -k 2>>$error >> $log\n\n";
            print OUT $bootTrees;
            if ( system($bootTrees) ) { die "Error running $bootTrees.\n"; }
            my $bestTree
            = "raxmlHPC-PTHREADS -p 10 -T $thread -f b -m GTRGAMMAI -t $outdir/RAxML_bestTree.$name -s $outdir/$name\_all_alignment.fna -z $outdir/RAxML_bootstrap.$name\_b -w $outdir -n $name\_best 2>>$error >> $log\n\n";
            print OUT $bestTree;
            if ( system($bestTree) ) { die "Error running $bestTree.\n"; }

        return "Bootstrap complete";
        close OUT;
    }
    }

    if ( $tree == 3 || $tree == 4 ) {
        print OUT "Reconstructing phylogeny using IQ-tree after finding the best model\n";
        if ($bsignal == 1){
            print OUT "Also bootstraping IQ-Trees trees\n";
            print OUT "\n";
            my $iqtree
                = "iqtree -m TEST -b $bootstrap -s $outdir/$name\_all_alignment.fna -nt $thread 2>>$error >> $log\n\n";
            print OUT $iqtree;
            if ( system($iqtree) ) { die "Error running $iqtree.\n"; }
        } else {
        my $iqtree
            = "iqtree -m TEST -s $outdir/$name\_all_alignment.fna -nt $thread 2>>$error >> $log\n\n";
        print OUT $iqtree;
        if ( system($iqtree) ) { die "Error running $iqtree.\n"; }
        }
    }
    print OUT "Tree phylogeny complete.\n";
    close OUT;

    if ( $tree == 1 ) { return ("Fasttree phylogeny complete"); }
    if ( $tree == 2 ) { return ("RAxML phylogeny complete"); }
    if ( $tree == 3 ) { return ("IQ-TREE Phylogeny complete"); }
    if ( $tree == 4 ) { return ("All phylogeny complete"); }
}
#------------------------------------------------------------------------------------------------------------------------------------------------#
sub bootstrap {
    my $bindir    = shift;
    my $outdir    = shift;
    my $thread    = shift;
    my $tree      = shift;
    my $name      = shift;
    my $bootstrap = shift;
    my $error     = shift;
    my $log       = shift;

# This was being run for tree=1, which only makes FastTree
# I think its OK to comment this out
# if ($tree==1){
#    my $bootTrees="raxmlHPC-PTHREADS -p 10 -T $thread -m GTRGAMMAI -b 10000 -t $outdir/$name\.fasttree -s $outdir/$name\_snp_alignment.fna -w $outdir -N $bootstrap -n $name\_b -k 2>>$error >> $log\n\n";
#    print $bootTrees;
#    if (system ($bootTrees)){die "Error running $bootTrees.\n";}

#    my $bestTree="raxmlHPC-PTHREADS -p 10 -T $thread -f b -m GTRGAMMAI -t $outdir/$name\.fasttree -s $outdir/$name\_snp_alignment.fna -z $outdir/RAxML_bootstrap.$name\_b -w $outdir -n $name\_best 2>>$error >> $log\n\n";
#    print $bestTree;
#    if (system ($bestTree)){die "Error running $bestTree.\n";}
# }
    open( OUT, ">>$log" );
    if ( $tree == 2 || $tree == 4 ) {
        my $bootTrees
            = "raxmlHPC-PTHREADS -p 10 -T $thread -m GTRGAMMAI -b 10000 -t $outdir/RAxML_bestTree.$name -s $outdir/$name\_all_alignment.fna -w $outdir -N $bootstrap -n $name\_b -k 2>>$error >> $log\n\n";
        print OUT $bootTrees;
        if ( system($bootTrees) ) { die "Error running $bootTrees.\n"; }
        my $bestTree
            = "raxmlHPC-PTHREADS -p 10 -T $thread -f b -m GTRGAMMAI -t $outdir/RAxML_bestTree.$name -s $outdir/$name\_all_alignment.fna -z $outdir/RAxML_bootstrap.$name\_b -w $outdir -n $name\_best 2>>$error >> $log\n\n";
        print OUT $bestTree;
        if ( system($bestTree) ) { die "Error running $bestTree.\n"; }

    }
    # if ( $tree == 3 || $tree == 4 ) {
    #     print OUT "Reconstructing phylogeny using IQ-tree after finding the best model\n";
    #     print OUT "Also bootstraping IQ-Trees trees\n";
    #     print OUT "\n";
    #     my $iqtree
    #         = "iqtree -redo -m TEST -b $bootstrap -s $outdir/$name\_all_alignment.fna -nt $thread 2 >>$error >> $log\n\n";
    #     print OUT $iqtree;
    #     if ( system($iqtree) ) { die "Error running $iqtree.\n"; }
    # }

    return "Bootstrap complete";
    close OUT;
}

sub extractGenes {
    my $dir      = shift;
    my $stat     = shift;
    my $file     = shift;
    my $bindir   = shift;
    my $list     = shift;
    my $thread   = shift;
    my $gapfile  = shift;
    my $genefile = shift;
    my $error    = shift;
    my $log      = shift;
    my $genedir  = $dir . '/PSgenes';

    open( OUT, ">>$log" );
    print OUT "\n";
    my $extract
        = "extractGenes.pl -d $dir -t $thread -l $list -s $stat -f $file -p $gapfile -g $genefile 2>>$error >> $log\n\n";
    print OUT $extract;
    if ( system($extract) ) { 
        die "Error running $extract Check $log and $error for details"; }

    opendir( DIR, "$genedir" ) || die "$!";
OUTER: while ( my $files = readdir(DIR) ) {
        next if ( $files =~ /^..?$/ );
        if ( $files =~ /.+\.fna$/ ) {
            my $file = $genedir . '/' . $files;
            open( IN, "$file" ) || die "$!";
            while (<IN>) {
                if ( !/^>/ ) {
                    if ( !/^ATG/ ) {
                        # print "$file\n";
                        `rm $file`;
                        next OUTER;
                    }
                }
            }
            close IN;
        }
    }
    close DIR;
    close OUT;
    return ("Genes with SNPs are now in PSgenes Directory");
}

sub translateGenes {
    my $dir     = shift;
    my $bindir  = shift;
    my $thread  = shift;
    my $program = shift;
    my $gencode = shift;
    my $error   = shift;
    my $log     = shift;
    my $genedir = $dir . '/PSgenes';

    open( OUT, ">>$log" );
    print OUT "\n";
    my $translate
        = "parallel_run.pl -d $genedir -t $thread -m $program -g $gencode 2>>$error >> $log\n\n";
    print OUT $translate;
    if ( system($translate) ) { die "Error running $translate.\n"; }
    close OUT;
    return ("Gene translation complete");
}

sub alignGenes {
    my $dir     = shift;
    my $bindir  = shift;
    my $thread  = shift;
    my $program = shift;
    my $gencode = shift;
    my $error   = shift;
    my $log     = shift;
    my $genedir = $dir . '/PSgenes';

    open( OUT, ">>$log" );
    print OUT "\n";
    my $align
        = "parallel_run.pl -d $genedir -t $thread -m $program -g $gencode 2>>$error >> $log\n\n";
    print OUT $align;
    if ( system($align) ) { die "Error running $align.\n"; }
    close OUT;
    return ("MSA complete");
}

sub revTransGenes {
    my $dir     = shift;
    my $bindir  = shift;
    my $thread  = shift;
    my $program = shift;
    my $gencode = shift;
    my $error   = shift;
    my $log     = shift;
    my $genedir = $dir . '/PSgenes';

    open( OUT, ">>$log" );
    print OUT "\n";    
    my $revTrans
        = "parallel_run.pl -d $genedir -t $thread -m $program -g $gencode 2>>$error >> $log\n\n";
    print OUT $revTrans;
    if ( system($revTrans) ) { die "Error running $revTrans.\n"; }
    close OUT;
    return ("Reverse translation complete,codon MSA complete");
}

sub core {
    my $dir     = shift;
    my $bindir  = shift;
    my $output  = shift;
    my $error   = shift;
    my $log     = shift;
    my $genedir = $dir . '/PSgenes';
    
    open( OUT, ">>$log" );
    print OUT "\n";   
    my $core = "catAlign.pl $genedir $output 2>>$error >> $log\n\n";
    print OUT $core;
    if ( system($core) ) { die "Error running $core.\n"; }
    close OUT;
    return ("Core gene alignment complete");
}

sub paml {
    my $dir     = shift;
    my $bindir  = shift;
    my $tree    = shift;
    my $model   = shift;
    my $suffix  = shift;
    my $NSsites = shift;
    my $core    = shift;
    my $thread  = shift;
    my $error   = shift;
    my $log     = shift;
    my $branch  = 1;
    my $pamldir = $dir . '/paml';

    open( OUT, ">>$log" );

    if ( $model == 0 ) {
        print OUT "\n";
        my $ps
            = "runPAML.pl -i $dir -t $thread -r $tree -m $model -n $NSsites -s $suffix 2>>$error >> $log\n\n";
        print OUT $ps;
        if ( system($ps) ) { die "Error running $ps.\n"; }
        my @site_files = glob("$pamldir/*/*$suffix");
        if (@site_files) { # @site_files is not empty...
            `mv $pamldir/*/*$suffix $pamldir`;
            print OUT "\n";
        } 

        my $parse
            = "parseSitePAML.pl $pamldir $NSsites 2>>$error >> $log\n\n";
        print OUT $parse;
        if ( system($parse) ) { die "Error running $parse. \n"; }
    }

    if ( $model == 2 ) {
        print OUT "\n";
        my $edit = "ParseTree.pl $tree 2>>$error >> $log\n\n";
        print OUT $edit;
        if ( system($edit) ) { die "Error running $edit.\n"; }

        # make the tree with branch label
        my($tree_name, $dirs, $tree_suffix) = File::Basename::fileparse($tree, qr/\.[^.]*/);
        my $branch_labeled_tree = $pamldir . "/$tree_name"."_BranchNumber" . "$tree_suffix";
        print OUT "\n";
        my $ps
            = "runPAML.pl -i $dir -t $thread -r $branch_labeled_tree -m $model -n $NSsites -s $suffix 2>>$error >> $log\n\n";
        print OUT $ps;
        if ( system($ps) ) { die "Error running $ps.\n"; }
        my @site_files = glob("$pamldir/*/*$suffix");
        if (@site_files) { # @site_files is not empty...
            `mv $pamldir/*/*$suffix $pamldir`;
            print OUT "\n";
        }

        my $parse
            = "parseSitePAML.pl $pamldir $NSsites 2>>$error >> $log\n\n";
        print OUT $parse;
        if ( system($parse) ) { die "Error running $parse. \n"; }

    }
    close OUT;
    return ("PAML run complete.");
}


sub PickRefGenome {

    # A function that picks the best genome to use as reference based on
    # ANI calculated from comparing sketches using mash (as implemented in
    # bbmap)
    my $workdir       = shift;
    my $refdir        = shift;
    my $error         = shift;
    my $log           = shift;
    my $sketch_output = shift;
    my $sketch_dir;
    my $ref_file;
    my %refs;
    my @references;
    my @contigs;
    my $contigfiles;

    open( OUT, ">>$log" );
    if ( $workdir =~ /.+\/$/ ) { my $temp = chop($workdir); }
    $sketch_dir = $workdir . '/sketches';
    mkdir $sketch_dir unless -d $sketch_dir;

    opendir( WORKDIR, $workdir ) or die $!;
    my @contigfiles
        = grep { /\.contig$|\.contigs$|\.ctg$]$/ && -f "$workdir/$_" }
        readdir(WORKDIR);
    my @fullcontig = map { $workdir . '/' . $_ } @contigfiles;
    closedir(WORKDIR);

    opendir( WORKDIR, $workdir ) or die $!;
    my @readfiles
        = grep { /\.fastq$|\.fq$]$/ && -f "$workdir/$_" } readdir(WORKDIR);

    my @forwardreads = grep !/R2\.fastq/, @readfiles;
    my @fullread = map { $workdir . '/' . $_ } @forwardreads;
    closedir(WORKDIR);

    opendir( REFDIR, $refdir ) or die $!;
    my @reffiles
        = grep { /\.fna$|\.fasta$|\.fa$|\.fsa$/ && -f "$refdir/$_" }
        readdir(REFDIR);
    my @fullref = map { $refdir . '/' . $_ } @reffiles;
    closedir(REFDIR);

    my @inputs = ( @fullcontig, @fullread, @fullref );

    # create sketch for all inputs
    foreach (@inputs) {
        my ( $filename, $path, $suffix )
            = File::Basename::fileparse( $_, qr/\.[^.]*/ );

        my $out_sketch = $sketch_dir . "/$filename" . ".sketch";
        print OUT "Sketching $_ \n";
        system(
            "sketch.sh in=$_ out=$out_sketch  mode=single name0=$_ fname=$_ name=$_ 2>>$error >> $log\n\n"
        );
    }

    # compare each sketch to each other
    opendir( SKDIR, $sketch_dir ) or die $!;
    my @sketchfiles
        = grep { /\.sketch$/ && -f "$sketch_dir/$_" } readdir(SKDIR);
    my @fullsketch = map { $sketch_dir . '/' . $_ } @sketchfiles;
    closedir(SKDIR);

    my $all_sketch = join ",", @fullsketch;

    my $ref_sketch = $sketch_dir . "/*.sketch";
    print OUT "Comparing sketches....\n";
    system(
        "comparesketch.sh alltoall compareself=f format=3 $sketch_dir/*.sketch| uniq | grep ^$refdir | sed 's/\t\t/\t/g'> $sketch_output"

    );

    my $ref_genome = `
    awk 'NR>1{
              arr[\$1] += \$3
              count[\$1] += 1
   }
   END{
      for (a in arr) {
         print a ',' arr[a] / count[a]}}' $sketch_output | sort -r -k 2 |head -n 1`;

   my $ref_genome1 = ( split / /, $ref_genome )[0];

    return $ref_genome1;
    close OUT;
}


# sub FilterGenomes {

#     # A function that outputs list of genomes that needs to be filtered out based on cutoff.
#     my $sketch_file   = shift;
#     my $ref_genome    = shift;
#     my $cutoff = shift;
#     my @remove_list;
    
#     $ref_genome =~ s/\//\\\//g;
#     my $x = `cat $sketch_file | grep $ref_genome |awk '{ if (\$3 < $cutoff ) print \$1, \$2 }'| sed 's/$ref_genome//g' | sed 's/ //g'`;
#     @remove_list = split(/\n/, $x);
#     return @remove_list;
# }


sub combo {
    my $by = shift;
    return sub { () }
        if !$by || $by =~ /\D/ || @_ < $by;
    my @list = @_;

    my @position = ( 0 .. $by - 2, $by - 2 );
    my @stop     = @list - $by .. $#list;
    my $end_pos  = $#position;
    my $done     = undef;

    return sub {
        return () if $done;
        my $cur = $end_pos;
        {
            if ( ++$position[$cur] > $stop[$cur] ) {
                $position[ --$cur ]++;
                redo if $position[$cur] > $stop[$cur];
                my $new_pos = $position[$cur];
                @position[ $cur .. $end_pos ] = $new_pos .. $new_pos + $by;
            }
        }
        $done = 1 if $position[0] == $stop[0];
        return @list[@position];
        }
}


sub SNPsAnalysis {
    # A function that parses .snps and .vcf file
    # it outputs another file with detail information on SNP position
    # like did SNP result in synonymous or non-synonymous changes
    # it requires a modified gff file that is generated from codingRegions()
    # the function runs chienchi's script that is implemented in EDGE
    # it only parses snps and vcf file of contigs and reads
    my $cds_gff = shift;
    my $snp_dir = shift;
    my $ref_fasta = shift;
    my $log = shift;
    my $error = shift;
    my $snp_file;

    open( OUT, ">>$log" );
    opendir( DIR, $snp_dir );
    while ( my $files = readdir(DIR) ) {
        next if ( $files =~ /^..?$/ );
        if ( $files =~ /contigs?.snps$/ ) {
            $snp_file = $snp_dir . '/' . $files;
            my $p = "SNP_analysis.pl -gff $cds_gff -SNP $snp_file -format nucmer -fasta $ref_fasta -output $snp_dir 2>>$error >> $log\n\n";
            print OUT $p;
            if ( system($p) ) { die "Error running $p.\n"; }
        }
        elsif ( $files =~ /.+\.vcf/ ) {
            $snp_file = $snp_dir . '/' . $files;
            my $p = "SNP_analysis.pl -gff $cds_gff -SNP $snp_file -format vcf -fasta $ref_fasta 2>>$error -output $snp_dir >> $log\n\n";
            print OUT $p;
            if ( system($p) ) { die "Error running $p.\n"; }
            }
    }
    close OUT;
}


# sub filter_genomes {

#     # A function that gives list of genomes whose ANI when compared to reference is less than the threshold.
#     # It needs a three column sketch file as input
#     # sketch file contains file path in specified ref directory
#     # but it outputs the path to folder where genomes are copied
#     my $sketch_file   = shift;
#     my $ref_genome    = shift;
#     my $cutoff = shift;
#     my $workdir = shift;
#     my @remove_list;
#     my @final_remove_list;
    
#     $ref_genome =~ s/\//\\\//g;
#     my $x = `cat $sketch_file | grep $ref_genome |awk '{ if (\$3 < $cutoff ) print \$1, \$2 }'| sed 's/$ref_genome//g' | sed 's/ //g'`;
#     @remove_list = split(/\n/, $x);

#     foreach my $remove_file (@remove_list){
#          my($filename, $dirs, $suffix) = File::Basename::fileparse($remove_file, qr/\.[^.]*/);
#          $filename =~ s/\W/_/g;
#          $filename = $filename . ".fna";
#          push(@final_remove_list, $filename)

#     }
#     return uniq(@final_remove_list);
# }

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


sub hyphy {
    my $dir      = shift;
    my $bindir   = shift;
    my $tree     = shift;
    my $outtree  = shift;
    my $core     = shift;
    my $threads  = shift;
    my $analysis = shift;
    my $error    = shift;
    my $log      = shift;

    $threads = int( $threads / 2 );

    my $hyphy
        = "runHyPhy.pl -i $dir -t $threads -r $tree -o $outtree -c $core -a bsrel 2>>$error >> $log\n\n";
    print $hyphy;
    if ( system($hyphy) ) { die "Error running $hyphy. \n"; }
}

1;
