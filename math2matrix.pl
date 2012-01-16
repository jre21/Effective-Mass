#!/usr/bin/perl

# Read from standard input an array in Mathematica's CForm, and output
# (to standard output) C code to fill an ARRAY (defined in defs.h)
# with its elements.
use strict;
use warnings;

my $line = "";
my @output = (0,0,"","");
my $i = 0;
my $j = 0;

while (<STDIN>)
{
    # Add next line of input to line, first removing whitespace and
    # trailing backslashes.
    s/\s*|\\$//g;
    $line .= $_;
}

# Convert Mathematica'isms to proper C code (some from GSL).
$line =~ s/Power/pow/g;
$line =~ s/Sqrt/sqrt/g;
$line =~ s/\.([^\d])/.0$1/g;
$line =~ s/Erfc/gsl_sf_erfc/g;
$line =~ s/Erf/gsl_sf_erf/g;
$line =~ s/Pi/M_PI/g;
$line =~ s/\bE\b/M_E/g;

while($line ne "")
{

    # Strip leading List( directives.  Unless this is the beginning of
    # input, increment row, reset column, and remove trailing
    # close-paren from output.
    while($line =~ /^List\((.*)$/)
    {
	$line = $1;
	if($output[2] ne "")
	{
	    if($output[2] =~ /(.*)\)/) { $output[2] = $1; }
	    if($output[3] =~ /(.*)\)/) { $output[3] = $1; }
	    $i++;
	    $j=0;
	}
    }

    # Format and print output.
    if($output[2] ne "")
    {
	print "gsl_matrix_complex_set(output, $output[0], $output[1], "
	    . "gsl_complex_rect($output[2], $output[3]));\n";
    }

    # Reset output.
    @output = ($i,$j,"0","0");

    # Continue parsing until we run out of commas, at which point we
    # have one remaining entry.
    if(!($line =~ /^[^,]*((Complex\([^\)]*\)|pow(\(((?:(?>[^()]+)|(?3))*)\)))[^,]*)*$/))
    {

	# If Complex is used, separate the real and imaginary parts, and
	# consume everything up to and including the complex number.
	if($line =~ /^([^,]*)Complex\(([^,\)]*),([^\)]*)\)(.*)$/)
	{
	    $output[2] = $1.$2;
	    $output[3] = $1.$3;
	    $line = $4;
	}
	# Otherwise blank the real part so it isn't assumed to be 0.
	else { $output[2] = ""; }

	# Check whether either output is 0 so we can avoid further
	# processing on it.
	if($output[2] =~ /^\(*0$/) { $output[2] = "0"; }
	if($output[3] =~ /^\(*0$/) { $output[3] = "0"; }

	# If Power is used, consume and add to output so we don't
	# parse the internal comma.
	while($line =~ /^([^,]*pow(\(((?:(?>[^()]+)|(?2))*)\)))(?'rem'.*)$/)
	{
	    if($output[2] ne "0") { $output[2] .= $1; }
	    if($output[3] ne "0") { $output[3] .= $1; }
	    $line = $+{rem};
	}

	# Add everything up to first comma to output and consume from
	# input.
	if($line =~ /^([^,]*),(.*)$/)
	{
	    $line = $2;
	    if($output[2] ne "0") { $output[2] .= $1; }
	    if($output[3] ne "0") { $output[3] .= $1; }
	}

	# Increment j.
	$j++;
    }

    # Finally, parse the last element.  This must be real due to
    # hermiticity, but we need to clear out two closing parentheses.
    else
    {
	if($line =~ /(.*)\)\)$/)
	{
	    print "gsl_matrix_complex_set(output, $i, $j, "
		. "gsl_complex_rect($1, 0));\n";
	}
	$line = "";
    }
}
