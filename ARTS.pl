#!/usr/bin/perl
# ARTS: Automated Randomization of multiple Traits for Study Design, using diploidly GA
# Mark Maienschein-Cline, last updated 8/19/2013
# mmaiensc@uic.edu
# Center for Research Informatics, University of Illinois at Chicago
#
# Copyright (C) 2013 Mark Maienschein-Cline
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.





use Getopt::Long qw(:config no_ignore_case);
#use Time::HiRes qw( clock_gettime );
use Math::Trig;
$|++;

#
# initialize random number parameters
#
&ran1_init();

#
# read command line
#
&read_command_line();

#
# read phenotype list: print the title lines of columns used for verbose
#
&read_data();
if( $verb eq "y" || $verb eq "l" ){
	printf("Using traits:");
	for($i=0; $i<= $#allcols; $i++){
		print "\t$titlevals[$allcols[$i]]";
	}
	print "\n";
	printf("Using trait combinations:");
	for($i=0; $i<= $#cols; $i++){
		printf("\t{%s", $titlevals[$cols[$i][0]]);
		for($j=1; $j<= $#{$cols[$i]}; $j++){
			printf(",%s", $titlevals[$cols[$i][$j]]);
		}
		printf("}");
	}
	print "\n";
}

#
# initialize GA parameters
#
&ga_init();

#
# if using the batchcolumn, fill in the batch
#
if( $bcolumn ne "" ){
	if( $verb eq "y" ){
		printf("Looking at column %i (%s) for batch assignments\n", $bcolumn+1, $titlevals[$bcolumn]);
	}
	# fill in batch from last column of @data
	$numbatches = 0;
	$foundbatchhash = {};
	@batchsizes = ();
	for($i=0; $i<=$#data; $i++){
		if( $foundbatchhash->{$data[$i][$bcolumn]} eq "" ){
			$foundbatchhash->{$data[$i][$bcolumn]} = $numbatches;
			$numbatches++;
			push(@batchnames, $data[$i][$bcolumn]);
			push(@batchsizes, 0);
		}
		$batchsizes[$foundbatchhash->{$data[$i][$bcolumn]}]++;
		$data[$i][$#{$data[0]}] = $data[$i][$bcolumn];
	}
	$mi = &mutual_info( $numbatches );
	$bestmi = $mi;
}

#
# else do sampling: run GA
#
if( $bcolumn eq "" ){
	&initialize_population();

	$oldavg = 1;
	$err = 0.0001;
	for($n=0; $n< $numgen; $n++){
		&add_immigrants();
		@population = &permute( \@population );
		$k = 0;
		$k = &crossover( $k );
		$k = &mutate( $k );
		$k = &add_parents( $k );
		@pool = sort{$a->{score} <=> $b->{score}} @pool;
		$average = &fill_population();

		# check if we've done enough already, and print out status
		if( $verb eq "y" ){printf("  Generation %i of %i, average fitness %0.4f\n", $n+1, $numgen, $average );}
		if( $oldavg >= $average && $oldavg - $average < $err ){last;}
		$oldavg = $average;
	}

	# save the final best one
	for($i=0; $i<= $#data; $i++){
		&fill_assignments( \@{$population[0]->{assignments}} );
	}
}

#
# print final log to stdout
#
if( $verb eq "y" || $verb eq "l" ){&print_info;}

#
# print result
#
if( $out ne "" ){
	open(OUT,">$out");
	print OUT "$title\t$bname\n";
	for($i=0; $i<= $#data; $i++){
		$orig[$i][1] = $data[$i][$#{$data[0]}];
		printf OUT ("%s\t%i\n", $orig[$i][0], $orig[$i][1]);
	}
	close(OUT);
}

###############
# SUBROUTINES #
###############

# read command line options
sub read_command_line{
	my $i;
	
	#
	# option default values
	#
	$in = "";
	$out = "";
	$bcolumn = "";
	$batch = "";
	$bname = "batch";
	$phenocols = "";
	$contcols = "";
	$datecols = "";
	$bins = 5;
	@blist = ();
	$verb = "y";
	$mmi = 0;

	$options = "
Usage: ./ARTS.pl <OPTIONS>
	REQUIRED:
	-i  input traits (rectangular, tab-delimited matrix, including title line with column names)
	-c  trait columns to randomize
	    comma- and semicolon delimited list, columns indexed from 1
	    all traits indicated by commas are used in joint distributions
	AND EITHER -b AND -o, OR -p:
	-b  batch sizes (a single number, or a comma-delimited list)
	-o  output file (formatted same as input file, with batch added as last column)
	   -or-
	-p  <batch column index>: print MI statistic for input traits using this column as batch designations
	    -p will not do any sampling
	OTHER OPTIONS:
	-cc continuously-valued columns (will be binned)
	-cd date-valued columns (should be M/D/Y); should also list these as continuous (in -cc)
	-cb number of bins to use for continuous or date columns (default: $bins for each)
	    can give 1 value, or a list of the same length as -cc; if a list, will be assigned in the same order as -cc
	-bn batch name (title of added column, default $bname)
	-s  random number seed (large negative integer, default: $seed)
";

#
# Secret options:
# -v y or l (verbose: print all, or just print status from beginning or end)
# -mmi force use of MMI objective function on all columns indicated by -c, over-riding any other settings from -c
#

	GetOptions('i=s'  => \$in,
		   'o=s'  => \$out,
		   'p=i'  => \$bcolumn,
		   'b=s'  => \$batch,
		   'c=s'  => \$phenocols,
		   'cc=s' => \$contcols,
		   'cd=s' => \$datecols,
		   'cb=s' => \$bins,
		   'bn=s' => \$bname,
		   's=i'  => \$seed,
		   'mmi'  => \$mmi,
		   'v=s'  => \$verb,
		) || die "$options\n";

	#
	# check that required inputs exist
	#
	if( $in eq "" ){&exit_required("i");}
	if( ($out eq "" || $batch eq "") && $bcolumn eq "" ){&exit_required("b and -o, or -p,");}
	if( $phenocols eq "" || $phenocols eq "None" ){&exit_required("c");}

	#
	# check that inputs values are OK
	#
	if( $bcolumn ne "" ){
		if( $bcolumn < 1 ){&exit_err("p","at least 1");}
		$bcolumn--;
	}
	if( $verb ne "y" && $verb ne "n" && $verb ne "l" ){&exit_err("v","y or n or l");}
	if( $seed > 0 ){$seed*= -1;}
	if( $seed == 0 ){&exit_err("s","non-zero");}

	#
	# if mmi, reset phenocols value using all found columns
	#
	if( $mmi ){
		@initcs = split(/[,;]/, $phenocols);
		# remove duplicates
		@clist = ();
		$cinds = {};
		for($i=0; $i<= $#initcs; $i++){
			if( $cinds->{$initcs[$i]} eq "" ){
				$cinds->{$initcs[$i]} = 1;
				push(@clist, $initcs[$i]);
			}
		}
		# add to new phenocols
		$phenocols = "$clist[0]";
		for($i=1; $i<= $#clist; $i++){
			$phenocols = sprintf("%s,%s", $phenocols, $clist[$i]);
		}
		$phenocols = sprintf("%s;%s", $phenocols, $clist[0]);
		for($i=1; $i<= $#clist; $i++){
			$phenocols = sprintf("%s;%s", $phenocols, $clist[$i]);
		}
	}

	#
	# extract phenotype columns
	#
	@cols = ();
	@allcols = ();
	$alllist = {};
	@jointlist = split(';',$phenocols);
	for($i=0; $i<= $#jointlist; $i++){
		@tmp = split(',',$jointlist[$i]);
		@tmp = &fix_cols( \@tmp );
		push(@cols, [@tmp]);
		for($j=0; $j<= $#tmp; $j++){
			if( $alllist->{$tmp[$j]} eq "" ){
				$alllist->{$tmp[$j]} = 1;
				push(@allcols, $tmp[$j]);
			}
		}
	}

	#
	# extract continuous and date columns
	# sort continuous columns so that bins correspond to them in order
	#
	if( $contcols ne "" && $contcols ne "None" ){
		@conts = split(',',$contcols);
		@conts = &fix_cols( \@conts );
		$numconts = $#conts+1;
	}
	if( $datecols ne "" && $datecols ne "None" ){
		@dates = split(',',$datecols);
		@dates = &fix_cols( \@dates );
		$numdates = $#dates+1;
		# check that date columns are among continuous columns
		for($i=0; $i<= $#dates; $i++){
			for($j=0; $j<= $#conts; $j++){
				if( $dates[$i] == $conts[$j] ){last;}
				if( $j==$#conts ){
					printf("Error: please specify date column %i as continuous\n", $dates[$i]+1 );
					die;
				}
			}
		}
	}
	if( $bins =~ /,/ ){
		@blist = split(',',$bins);
		if( $#blist+1 != $#conts + 1 ){
			printf("Error: you input %i bins, but %i columns that need binning\n", $#blist+1, $#conts+1);
			die;
		}
	}
	else{
		for($i=0; $i<= $#conts; $i++){
			push(@blist, $bins);
		}
	}		
}

# print error message for flag $_[0], with correct values $_[1], and print usage
sub exit_err{
        printf("Error: set -%s to be %s\n%s\n", $_[0], $_[1], $options);
        exit;
}
# print error message saying flag $_[0] is required
sub exit_required{
        printf("Error: option -%s is required\n%s\n", $_[0], $options);
        exit;
}

# fix all indices in array $_[0]: cast to integer, check at least 1, and subtract 1
sub fix_cols{
	my @list;
	my $i;
	@list = @{$_[0]};
	for($i=0; $i<= $#list; $i++){
		$list[$i] = sprintf("%i", $list[$i]);
		if( $list[$i] < 1 ){
			print "Error: column indices should be at least 1\n";
			die;
		}
		$list[$i]--;
	}
	return @list;
}

# print info about best randomization
sub print_info{
	#
	# get MI of each phenotype and average
	#
	$bestmi = &mutual_info();
	@bestmilist = &individual_mi( $numbatches );
	$bestavgmi = 0;
	for($i=0; $i<= $#bestmilist; $i++){
		$bestavgmi+= $bestmilist[$i]/($#bestmilist+1);
	}

	printf("Final MI %0.4f ; Individual trait MIs (mean %0.4f ): ", $bestmi, $bestavgmi);
	for($i=0; $i<= $#bestmilist; $i++){
		printf("\t%0.4f", $bestmilist[$i]);
	}
	print "\n-----------------------------------------------------------------\n";
	# 
	# print the counts for each phenotype in each batch
	# 
	# first title line: phenotype names
	for($i=0; $i<= $#allcols; $i++){
		printf("\t%s values", $titlevals[$allcols[$i]]);
		for($j=1; $j<= $#{$items->{$allcols[$i]}->{list}}; $j++){
			printf("\t");
		}
	}
	print "\nBatch (size)";
	# second title line: phenotype values
	for($i=0; $i<= $#allcols; $i++){
		for($j=0; $j<= $#{$items->{$allcols[$i]}->{list}}; $j++){
			if( $items->{$allcols[$i]}->{list}[$j] ne "" ){printf("\t%s", &name($items->{$allcols[$i]}->{list}[$j], $allcols[$i]) );}
			else{printf("\tempty");}
		}
	}
	print "\n-------";
	# print a line of dashes to separate
	for($i=0; $i<= $#allcols; $i++){
		for($j=0; $j<= $#{$items->{$allcols[$i]}->{list}}; $j++){
			printf("\t-------");
		}
	}
	print "\n";
	for($k=0; $k< $numbatches; $k++){
		printf("%s (%i)", $batchnames[$k], $batchsizes[$k] );
		for($i=0; $i<= $#allcols; $i++){
			for($j=0; $j<= $#{$items->{$allcols[$i]}->{list}}; $j++){
				printf("\t%i", &count( $#{$data[0]}, $batchnames[$k], $allcols[$i], $items->{$allcols[$i]}->{list}[$j] ) );
			}
		}
		print "\n";
	}
	print "-------";
	# print a line of dashes to separate
	for($i=0; $i<= $#allcols; $i++){
		for($j=0; $j<= $#{$items->{$allcols[$i]}->{list}}; $j++){
			printf("\t-------");
		}
	}
	# print totals for each type
	print "\nTotal";
	for($i=0; $i<= $#allcols; $i++){
		for($j=0; $j<= $#{$items->{$allcols[$i]}->{list}}; $j++){
			printf("\t%i", $items->{$allcols[$i]}->{$items->{$allcols[$i]}->{list}[$j]}[1] );
		}
	}
	print "\n";
}

# for continuous valued columns, checked by $_[1], convert value $_[0] back to a range
sub name{
	my $i;
	my $binw;

	# 
	# if there aren't any continuous columns, or $_[1] doesn't match one, just return $_[0]
	#
	if( $#conts < 0 ){return $_[0];}
	for($i=0; $i<= $#conts; $i++){
		if( $_[1] == $conts[$i] ){last;}
		if( $i==$#conts ){return $_[0];}
	}

	#
	# convert bin value back to continuous value
	#
	$binw = ($contstats[$i][2]-$contstats[$i][0])/$blist[$i];
	$val1 = $binw*$_[0]+$contstats[$i][0];
	$val2 = $binw*($_[0]+1)+$contstats[$i][0];
	
	# 
	# if there aren't any date columns, or $_[1] doesn't match one, just return the range val1-val2
	#
	if( $#dates < 0 ){return sprintf("%s-%s", $val1, $val2);}
	for($i=0; $i<= $#dates; $i++){
		if( $_[1] == $dates[$i] ){last;}
		if( $i==$#dates ){return sprintf("%s-%s", $val1, $val2);}
	}
	$val1 = sprintf("%i", $val1);
	$val2 = sprintf("%i", $val2);
	return sprintf("%s-%s", &convert_date( $val1 ), &convert_date( $val2 ));

}

# read in regular matrix from $in
# for continuous (including date-value) columns, make histograms
sub read_data{
	my @lines;
	my $i;
	@data = ();
	$items = {};
	@orig = ();
	@titlevals = ();
	@batchsizes = ();
	$numbatches = 0;

	#
	# fix newline convention
	#
	
	
	open(IN,"$in") || die "Error: can't open $in\n";
	#
	# read in all lines and check formatting
	#
	@lines = <IN>;
	if( $lines[0] =~ /\r/ && $#lines == 0 ){
		# this happens with tab-delimited text saved from excel
		@lines = split('\r', $lines[0]);
	}

	#
	# read title line
	#
	$title = $lines[0];
	chomp($title);
	@titlevals = split('\t',$title);
	for($k=1; $k<= $#lines; $k++){
		$line = $lines[$k];
		chomp($line);
		if( $line ne "" ){	# ignore blank lines
			@parts = split('\t',$line);
			for($i=$#parts+1; $i<= $#titlevals; $i++){
				push(@parts, "");
			}
			if( $#parts != $#titlevals ){
				printf("Error: not enough columns in line %i\n", $#data+2);
				die;
			}
			# push 1 extra for the batch
			push(@parts, 0);
			push(@data, [(@parts)] );
			push(@orig, [($line, 0)] );
		}
	}
	close(IN);

	#
	# exit if no data read
	#
	if( $#data < 0 ){
		printf("Error: no samples were read in\n");
		die;
	}

	#
	# if batch is not empty, check batches:
	#   if no commas, cast to integer and count how many are needed
	#   if there are commas, get batched on given sizes
	#   double check that we add up
	# 
	@batchnames = ();
	if( $batch ne "" ){
		if( $batch !~ /,/ ){
			$batch = sprintf("%i", $batch);
			# fix batch size if too big
			if( $batch > $#data + 1 ){
				printf("Warning: you have %i samples, but asked for a batch size of %i, so there is only 1 batch\n", $#data+1, $batch);
				$batch = $#data+1;
			}
			$numbatches = ($#data+1)/$batch;
			$exactbatches = sprintf("%i", $numbatches);
			if( $exactbatches < $numbatches ){$exactbatches++;}
			$numbatches = $exactbatches;
			for($i=0; $i< $numbatches-1; $i++){
				push(@batchsizes, $batch);
			}
			push(@batchsizes, $batch - ($numbatches*$batch - ($#data+1)) );
		}
		else{
			@batchsizes = split(',',$batch);
			$numbatches = $#batchsizes+1;
		}
		$tot = 0;
		for($i=0; $i< $numbatches; $i++){
			push(@batchnames, $i+1);
			$tot+= $batchsizes[$i];
		}
		if( $tot != $#data+1 ){
			printf("Error: have %i spaces in all batches, but %i samples\n", $tot, $#data+1);
			die;
		}
	}

	#
	# convert dates to numbers
	#
	for($i=0; $i<= $#data; $i++){
		for($j=0; $j<= $#dates; $j++){
			if( $data[$i][$dates[$j]] ne "" ){$data[$i][$dates[$j]] = &convert_date( $data[$i][$dates[$j]] );}
		}
	}

	# 
	# for all continuous columns, compute median and fill in missing values
	# also record max and min for binning
	#
	@contstats = ();	# records min, median, max for each continuous column
	for($j=0; $j<= $#conts; $j++){
		@tmp = ();
		for($i=0; $i<= $#data; $i++){
			if( $data[$i][$conts[$j]] ne "" ){push(@tmp, $data[$i][$conts[$j]]);}
		}
		@tmp = sort{ $a <=> $b } @tmp;
		$median = $tmp[sprintf("%i", ($#tmp+1)/2)];
		push(@contstats, [($tmp[0], $median, $tmp[$#tmp])] );
	#	for($i=0; $i<= $#data; $i++){
	#		if( $data[$i][$conts[$j]] eq "" ){$data[$i][$conts[$j]] = $median;}
	#	}
	}

	#
	# for all continuous columns, bin data
	#
	for($j=0; $j<= $#conts; $j++){
		$binw = ($contstats[$j][2] - $contstats[$j][0])/($blist[$j]);
		if( $binw == 0 ){
			printf("Error: max and min of column %i are equal (max/min are %f/%f)\n", $conts[$j]+1, $contstats[$j][2], $contstats[$j][0] );
			die;
		}
		for($i=0; $i<= $#data; $i++){
			if( $data[$i][$conts[$j]] ne "" ){
				$data[$i][$conts[$j]] = sprintf("%i", ($data[$i][$conts[$j]] - $contstats[$j][0])/$binw);
				if( $data[$i][$conts[$j]] >= $blist[$j] ){$data[$i][$conts[$j]] = $blist[$j]-1;}
			}
		}
	}

	# 
	# for each column we're using, count how many item types there are
	# empty phenotypes are considered their own, distinct phenotype
	#
	$items = &itemize( \@allcols );
}

# count how many item types of @{$_[0]} there are in @data
sub itemize{
	my $i;
	my $j;
	my $info;
	my @cols;
	my $items;
	@cols = @{$_[0]};

	for($j=0; $j<= $#cols; $j++){
		$info = {};
		$info->{list} = ();
		for($i=0; $i<= $#data; $i++){
			if( $info->{$data[$i][$cols[$j]]} eq "" ){
				$info->{$data[$i][$cols[$j]]} = [($#{$info->{list}}+1,0)];
				push(@{$info->{list}}, $data[$i][$cols[$j]]);
			}
			$info->{$data[$i][$cols[$j]]}[1]++;
			$info->{count}++;
		}
		$info->{num} = $#{$info->{list}}+1;
		$items->{$cols[$j]} = $info;
	}

	for($j=0; $j<= $#cols; $j++){
		# this set prints the number of values and counts for each phenotype
		#printf("%i,%s:", $cols[$j], $titlevals[$cols[$j]]);
		#for($k=0; $k<= $#{$items->{$cols[$j]}->{list}}; $k++){
		#	printf("\t%s,%i", $items->{$cols[$j]}->{list}[$k], $items->{$cols[$j]}->{$items->{$cols[$j]}->{list}[$k]}[1]  );
		#}
		#print "\n";
		if( $items->{$cols[$j]}->{num} > 20 ){
			printf("Warning: column %i (%s) has %i values; should you make it continuous?\n", $cols[$j], $titlevals[$cols[$j]], $items->{$cols[$j]}->{num} );
		}
	}
	return $items;
}

# convert date in M/D/Y to integer, or integer to M/D/Y
sub convert_date{
	my $date;
	my $month;
	my $day;
	my $year;
	my $months;
	my $i;
	# cumulative days per month
	$months->{0} = 0;
	$months->{1} = 31;
	$months->{2} = 59;
	$months->{3} = 90;
	$months->{4} = 120;
	$months->{5} = 151;
	$months->{6} = 181;
	$months->{7} = 212;
	$months->{8} = 243;
	$months->{9} = 273;
	$months->{10} = 304;
	$months->{11} = 334;
	$months->{12} = 365;
	$date = $_[0];

	# convert date to integer
	if( $date =~ /\// ){
		($month, $day, $year) = split('/',$date);
		$month = sprintf("%i", $month);
		$day = sprintf("%i", $day);
		$year = sprintf("%i", $year);
		if( $month < 1 || $month > 12 ){
			print "Error: found a month $month not between 1 and 12\n";
			die;
		}
		if( $day < 1 || $day > 31 ){
			print "Error: found a day $day not between 1 and 31\n";
			die;
		}
			
		return $day + $months->{$month-1} + $year*$months->{12};
	}
	# convert integer to date
	elsif( $date == sprintf("%i", $date) ){
		$year = sprintf("%i", $date/($months->{12}));
		$month = $date-$year*$months->{12};
		for($i=1; $i<=12; $i++){
			if( $month < $months->{$i} ){last;}
		}
		$day = $month - $months->{$i-1};
		$month = $i;
		return sprintf("%s/%s/%s", $month, $day, $year);
	}
	else{
		printf("\nError: unrecognized format in convert_date(): %s\n", $date);
		die;
	}
}

# set globals used by ran1
sub ran1_init{
	#
	# random number variables
	#
	$iset = 0;
	$gset = 0;
	#$iseed = clock_gettime(CLOCK_REALTIME);
	#($first, $second) = split('\.', $iseed);
	#$seed = sprintf("-%i%i", $second, $first);
	$seed = -10854829;
	$M1 = 259200;
	$IA1 = 7141;
	$IC1 = 54773;
	$RM1 = 1.0/$M1;
	$M2 = 134456;
	$IA2 = 8121;
	$IC2 = 28411;
	$RM2 = 1.0/$M2;
	$M3 = 243000;
	$IA3 = 4561;
	$IC3 = 51349;
	$iff = 0;
	$ix1 = 0;
	$ix2 = 0;
	$ix3 = 0;
	@ranarray = ();
	for($i=0; $i< 98; $i++){
		push(@ranarray, 0);
	}
}

# uniform random number generator, seed, iff, and various capital-letter variables set in beginning
sub ran1{
	my $j;
	my $temp;

	if( $seed < 0 || $iff == 0 ){
		$iff = 1;
		$ix1 = ($IC1 - $seed)%$M1;
		$ix1 = ($IA1*$ix1 + $IC1)%$M1;
		$ix2 = $ix1%$M2;
		$ix1 = ($IA1*$ix1 + $IC1)%$M1;
		$ix3 = $ix1%$M3;
		for($j=1; $j<= 97; $j++){
			$ix1 = ($IA1*$ix1 + $IC1)%$M1;
			$ix2 = ($IA2*$ix2 + $IC2)%$M2;
			$ranarray[$j] = ($ix1 + $ix2*$RM2)*$RM1;
		}
		$seed = 1;
	}
	$ix1 = ($IA1*$ix1 + $IC1)%$M1;
	$ix2 = ($IA2*$ix2 + $IC2)%$M2;
	$ix3 = ($IA3*$ix3 + $IC3)%$M3;

	$j = sprintf("%i", 1 + ((97*$ix3)/$M3) );
	if( $j> 97 || $j< 1 ){
		printf("Error in ran1: $j outside of [1:97]\n");
		die;
	}
	$temp = $ranarray[$j];
	$ranarray[$j] = ($ix1 + $ix2*$RM2)*$RM1;
	return $temp;
}

# permute array $_[0]
sub permute{
	my @assignments;
	my $i;
	my $j;
	my $tmp;

	@assignments = @{$_[0]};

	# 
	# shuffle batches randomly
	#
	for($i=$#assignments; $i>= 0; $i--){
		$j = sprintf("%i", ($i+1)*&ran1() );
		$tmp = $assignments[$j];
		$assignments[$j] = $assignments[$i];
		$assignments[$i] = $tmp;
	}
	return @assignments;
}

# fill data with assignments $_[0]
sub fill_assignments{
	my @list;
	my $i;
	@list = @{$_[0]};
	if( $#list != $#data ){
		print "Error in fill_assignments: mismatching list lengths\n";
		die;
	}
	for($i=0; $i<= $#list; $i++){
		$data[$i][$#{$data[0]}] = $list[$i];
	}
}

# compute mutual information of a batch assignment
sub mutual_info{
	my $i;
	my $s;
	my $mi;
	my $stot;

	$mi = 0;
	for($i=0; $i<= $#cols; $i++){
		$mi += &this_mi( $_[0], $#{$data[0]}, \@{$cols[$i]} )/($#cols+1);
	}
	
	return $mi;
}

# compute all single-phenotype mutual information
sub individual_mi{
	my $i;
	my @list;
	my @milist;
	
	@milist = ();
	for($i=0; $i<= $#allcols; $i++){
		@list = ($allcols[$i]);
		push(@milist, &this_mi( $_[0], $#{$data[0]}, \@list ) );
	}
	return @milist;
}

# compute mutual information of columns $_[1] ($_[0] bins) and all of @{$_[2]}
sub this_mi{
	my $i;
	my $j;
	my $summand;
	my @list;
	my $jprob;
	my $m1prob;
	my $m2prob;
	my $jbin;
	my $m1bin;
	my $m2bin;
	my $jbinstot;
	my $m1binstot;
	my $m2binstot;
	my @jbinlist;
	my @m1binlist;
	my @m2binlist;
	my $mi;
	my $s1;
	my $s2;
	my $s;
	@list = @{$_[2]};
	
	# initialize probabilities
	$jprob = {};			# joint distribution
	$m1prob = {};			# batch marginal dist
	$m2prob = {};			# pheno marginal dist
	@jbinlist = ();			# phenotype combos found in joint distribution
	@m1binlist = ();		# batches found in batches distribution (1st marginal dist)
	@m2binlist = ();		# phenotype combos found in phenotypes distribution (2nd marginal dist)

	#	
	# read through data and add to distributions
	#
	$summand = 1.0/($#data+1);
	for($i=0; $i<= $#data; $i++){
		#
		# define bin names based on phenotype/batch
		# for phenotypes p1, p2, etc., batch b:
		#   joint = p1_p2_..._pn_b
		#   1st marginal = b
		#   2nd marginal = p1_p2_..._pn
		#
		$jbin = sprintf("%s", $data[$i][$#{$data[0]}]);
		$m1bin = sprintf("%s", $data[$i][$#{$data[0]}]);
		$m2bin = "";
		for($j=0; $j<= $#list; $j++){
			# NOTE:
			# $list[$j] is a phenotype column (e.g., gender)
			# $data[$i][$list[$j]] is the value of that phenotype in sample $i (e.g., M or F)
			# $items->{$list[$j]}->{$data[$i][$list[$j]]}[0] is the bin index (e.g., M->0, F->1) of that phenotype

			$jbin = sprintf("%s_%i", $jbin, $items->{$list[$j]}->{$data[$i][$list[$j]]}[0]);
			if( $j>0 ){$m2bin = sprintf("%s_", $m2bin);}
			$m2bin = sprintf("%s%i", $m2bin, $items->{$list[$j]}->{$data[$i][$list[$j]]}[0]);
		}

		# 
		# check if we've already seen this bin, for each distribution
		# initialize probabilities and add to list if it's the first time
		#
		if( $jprob->{$jbin} eq "" ){
			$jprob->{$jbin} = 0;
			push(@jbinlist, [($jbin, $m1bin, $m2bin)] );
		}
		if( $m1prob->{$m1bin} eq "" ){
			$m1prob->{$m1bin} = 0;
			push(@m1binlist, $m1bin);
		}
		if( $m2prob->{$m2bin} eq "" ){
			$m2prob->{$m2bin} = 0;
			push(@m2binlist, $m2bin);
		}

		#
		# add a count to each distribution
		#
		$jprob->{$jbin} += $summand;
		$m1prob->{$m1bin} += $summand;
		$m2prob->{$m2bin} += $summand;
	}

	#
	# compute mutual information, and entropy of m1prob and m2prob (for normalization)
	#
	$mi = 0;
	$s1 = 0;
	$s2 = 0;
	for($i=0; $i<= $#jbinlist; $i++){
		$mi+= ($jprob->{$jbinlist[$i][0]}) * log( ($jprob->{$jbinlist[$i][0]})/($m1prob->{$jbinlist[$i][1]} * $m2prob->{$jbinlist[$i][2]}) );
	}
	for($i=0; $i<= $#m1binlist; $i++){
		$s1-= $m1prob->{$m1binlist[$i]} * log( $m1prob->{$m1binlist[$i]} );
	}
	for($i=0; $i<= $#m2binlist; $i++){
		$s2-= $m2prob->{$m2binlist[$i]} * log( $m2prob->{$m2binlist[$i]} );
	}
	$s = sqrt($s1*$s2);

	#	
	# normalize mi
	#
	if( $s>0 ){$mi/= $s;}

	#
	# return normalized mi (0=independent, 1=completely dependent)
	#
	return $mi;
}

# count how many of @data have column $_[0] equal $_[1] and column $_[2] equal $_[3]
sub count{
	my $i;
	my $tot;

	$tot = 0;
	for($i=0; $i<= $#data; $i++){
		if( $data[$i][$_[0]] eq $_[1] && $data[$i][$_[2]] eq $_[3] ){$tot++;}
	}
	return $tot;
}

# initialize GA parameters and large matrices
sub ga_init{
	my $i;
	my $j;
	my $info;

	$popsize = 100;
	$numgen = 300;
	$nchrmuts = 2;
	$nnewimm = 10;
	$nkeepparents = 2;
	$nchrpool = ($nnewimm+$popsize) + ($nnewimm+$popsize)*$nchrmuts + $nkeepparents;

	#
	# population to turn over each generation
	#
	@population = ();
	for($i=0; $i< $popsize+$nnewimm; $i++){
		$info = {};
		$info->{score} = 0;
		$info->{assignments} = [()];
		for($j=0; $j<= $#data; $j++){
			push(@{$info->{assignments}}, 0);
		}
		push(@population, $info);
	}

	#
	# new individuals to fill each generation
	#
	@pool = ();
	for($i=0; $i< $nchrpool; $i++){
		$info = {};
		$info->{score} = 0;
		$info->{assignments} = [()];
		for($j=0; $j<= $#data; $j++){
			push(@{$info->{assignments}}, 0);
		}
		push(@pool, $info);
	}

	#
	# array to randomize for batch assignments
	#
	@batched = ();
	@bcounts = ();
	for($i=0; $i<= $#batchsizes; $i++){
		push(@bcounts, 0);
		for($j=0; $j< $batchsizes[$i]; $j++){
			push(@batched, $i+1);
		}
	}
}

# initialize the population array: randomize $popsize batches and score each one
sub initialize_population{
	my $i;
	
	for($i=0; $i< $popsize; $i++){
		@{$population[$i]->{assignments}} = &permute( \@batched );
		&fill_assignments( \@{$population[$i]->{assignments}} );
		$population[$i]->{score} = &mutual_info( $numbatches );
	}
}

# complete the crossover step, keeping track of our index with $_[0]
sub crossover{
	my $i;
	my $j;
	my $k;

	$k = $_[0];

	for($i=0; $i< $popsize + $nnewimm; $i+=2){
		&do_cross($i, $k);
		&fill_assignments( \@{$pool[$k]->{assignments}} );
		$pool[$k]->{score} = &mutual_info( $numbatches );
		$k++;
		&fill_assignments( \@{$pool[$k]->{assignments}} );
		$pool[$k]->{score} = &mutual_info( $numbatches );
		$k++;
	}
	return $k;
}

# do a crossover between population members $_[0] and $_[0]+1, fill to pool members $_[1] and $_[1]+1
sub do_cross{
	my $i;
	my $j;
	my $popmem;
	my $poolmem;
	my $index;
	my @swap;
	my @subswap;
	my $swapinds;

	$popmem = $_[0];
	$poolmem = $_[1];
	$index = sprintf("%i", $#data * &ran1() );

	#
	# count how many of each batch to switch
	#
	for($i=0; $i<= $#bcounts; $i++){
		$bcounts[$i] = 0;
	}
	for($i=0; $i<= $index; $i++){
		$bcounts[$population[$popmem]->{assignments}[$i] - 1]++;
	}

	#
	# for each batch i:
	# 	record into subswap all indices with that batch from popmem+1
	#	permute subswap
	#	add first bcounts[$i] into swap
	# then sort swap
	#
	@swap = ();
	$swapinds = {};
	for($i=0; $i<= $#bcounts; $i++){
		@subswap = ();
		for($j=0; $j<= $#data; $j++){
			if( $population[$popmem+1]->{assignments}[$j] == $i+1 ){
				push(@subswap, $j);
			}
		}
		@subswap = &permute( \@subswap );
		for($j=0; $j< $bcounts[$i]; $j++){
			push(@swap, $subswap[$j]);
			$swapinds->{$subswap[$j]} = 1;
		}
	}
	@swap = sort{$a <=> $b} @swap;

	#
	# fill start of first new chr from swap indices of second chr, and end from end of first chr
	#
	for($i=0; $i<= $index; $i++){
		$pool[$poolmem]->{assignments}[$i] = $population[$popmem+1]->{assignments}[$swap[$i]];
	}
	for($i=$index+1; $i<= $#data; $i++){
		$pool[$poolmem]->{assignments}[$i] = $population[$popmem]->{assignments}[$i];
	}

	#
	# fill start of second chr from start of first chr, and end from remaining parts of second chr
	#
	for($i=0; $i<= $index; $i++){
		$pool[$poolmem+1]->{assignments}[$i] = $population[$popmem]->{assignments}[$i];
	}
	$j = $index+1;
	for($i=0; $i<= $#data; $i++){
		if( $swapinds->{$i} != 1 ){
			$pool[$poolmem+1]->{assignments}[$j] = $population[$popmem+1]->{assignments}[$i];
			$j++;
		}
	}
	

	#
	# check that batch counts are still ok
	#
	$checkbatch = 0;
	if( $checkbatch ){
		for($i=0; $i<= $#bcounts; $i++){
			$bcounts[$i] = 0;
		}
		for($i=0; $i<= $#data; $i++){
			$bcounts[$pool[$poolmem]->{assignments}[$i] - 1]++;
		}
		for($i=0; $i<= $#bcounts; $i++){
			if( $bcounts[$i] != $batchsizes[$i] ){
					print "Error in do_cross: lost some batch counts in first daughter chr\n";
				exit;
			}
		}
		for($i=0; $i<= $#bcounts; $i++){
			$bcounts[$i] = 0;
		}
		for($i=0; $i<= $#data; $i++){
			$bcounts[$pool[$poolmem+1]->{assignments}[$i] - 1]++;
		}
		for($i=0; $i<= $#bcounts; $i++){
			if( $bcounts[$i] != $batchsizes[$i] ){
				print "Error in do_cross: lost some batch counts in first daughter chr\n";
				exit;
			}
		}
	}
}

# complete the mutation step, keeping track of our index with $_[0]
sub mutate{
	my $i;
	my $j;
	my $k;

	$k = $_[0];

	for($i=0; $i< $popsize+$nnewimm; $i++){
		for($j=0; $j< $nchrmuts; $j++){
			&do_mutation($i, $k);
			&fill_assignments( \@{$pool[$k]->{assignments}} );
			$pool[$k]->{score} = &mutual_info( $numbatches );
			$k++;
		}
	}
	return $k;
}

# do a mutation for population member $_[0], fill to pool member $_[1]
sub do_mutation{
	my $i;
	my $popmem;
	my $poolmem;
	my $index1;
	my $index2;

	$popmem = $_[0];
	$poolmem = $_[1];

	#
	# fill all of poolmem
	#
	for($i=0; $i<= $#data; $i++){
		$pool[$poolmem]->{assignments}[$i] = $population[$popmem]->{assignments}[$i];
	}

	# 
	# switch two members
	#
	$index1 = sprintf("%i", ($#data+1) * &ran1() );
	$index2 = sprintf("%i", ($#data+1) * &ran1() );

	$pool[$poolmem]->{assignments}[$index1] = $population[$popmem]->{assignments}[$index2];
	$pool[$poolmem]->{assignments}[$index2] = $population[$popmem]->{assignments}[$index1];
}

# add immigrants, keeping track of our index with $_[0]
sub add_immigrants{
	my $j;

	for($j=0; $j< $nnewimm; $j++){
		@{$population[$popsize+$j]->{assignments}} = &permute( \@batched );
		&fill_assignments( \@{$population[$popsize+$j]->{assignments}} );
		$population[$popsize+$j]->{score} = &mutual_info( $numbatches );
		$k++;
	}
}

# add top-scoring parents, keeping track of our index with $_[0]
sub add_parents{
	my $i;
	my $j;
	my $k;

	$k = $_[0];
	# sort population now
	@population = sort{$a->{score} <=> $b->{score}} @population;

	for($j=0; $j< $nkeepparents; $j++){
		&copy_parents( $j, $k );
		# will copy score also in copy_parents()
		$k++;
	}
	return $k;
}

# add population member $_[0] to pool member $_[1]
sub copy_parents{
	my $i;
	my $popmem;
	my $poolmem;
	my $index1;
	my $index2;

	$popmem = $_[0];
	$poolmem = $_[1];

	#
	# fill all of poolmem
	#
	for($i=0; $i<= $#data; $i++){
		$pool[$poolmem]->{assignments}[$i] = $population[$popmem]->{assignments}[$i];
	}
	$pool[$poolmem]->{score} = $population[$popmem]->{score};
}

# copy top of pool to population (assumes pool is sorted)
sub fill_population{
	my $i;
	my $j;
	my $m;

	$m = 0;
	for($i=0; $i< $popsize; $i++){
		$population[$i]->{score} = $pool[$i]->{score};
		$m+= $population[$i]->{score};
		for($j=0; $j<= $#data; $j++){
			$population[$i]->{assignments}[$j] = $pool[$i]->{assignments}[$j];
		}
	}
	return $m/$popsize;
}




