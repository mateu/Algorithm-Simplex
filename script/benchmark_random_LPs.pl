#!/usr/bin/env perl
use strict;
use warnings;
use Benchmark;
use PDL::Lite;
use Getopt::Long;
use Algorithm::Simplex::Float;
use Algorithm::Simplex::PDL;
use Algorithm::Simplex::Rational;

=head1 Name

benchmark_random_LPs.pl - Benchmark the three models w/ random Linear Programs

=head1 Usage

perl benchmark_random_LPs.pl --rows 50 --columns 50 -n 50

=cut

my $rows          = 20;
my $columns       = 20;
my $number_of_LPs = 20;

GetOptions(
    'rows|r=i'          => \$rows,
    'columns|c=i'       => \$columns,
    'number_of_LPs|n=i' => \$number_of_LPs,
);

srand;
my $matrix = random_float_matrix( $rows, $columns, 1 );

# Get shell tableau object for access to EPSILON and MAXIMUM_PIVOTS
my $tableau_shell = Algorithm::Simplex->new( tableau => [ [] ] );

timethese(
    $number_of_LPs,
    {
        float    => 'solve_LP("float")',
        piddle   => 'solve_LP("piddle")',
        rational => 'solve_LP("rational")',
    }
);

=head2 solve_LP

Solver subroutine for a given model and initial tableau.

=cut

sub solve_LP {
    my $model   = shift;
    my $tableau = matrix_copy($matrix);

    my $tableau_object =
      $model eq 'float'
      ? Algorithm::Simplex::Float->new( tableau => $tableau )
      : $model eq 'piddle'
      ? Algorithm::Simplex::PDL->new( tableau => $tableau )
      : $model eq 'rational'
      ? Algorithm::Simplex::Rational->new( tableau => $tableau )
      : die "The model type: $model could not be found.";

    my $counter = 1;
    until ( $tableau_object->is_optimal ) {
        my ( $pivot_row_number, $pivot_column_number ) =
          $tableau_object->determine_bland_pivot_row_and_column_numbers;
        $tableau_object->pivot( $pivot_row_number, $pivot_column_number );
        $tableau_object->exchange_pivot_variables( $pivot_row_number,
            $pivot_column_number );
        $counter++;

        # Too many pivots?
        if ( $counter > $tableau_shell->MAXIMUM_PIVOTS ) {
            warn "HALT: Exceeded the maximum number of pivots allowed: "
              . $tableau_shell->MAXIMUM_PIVOTS . "\n";
            return 0;
        }
    }
    return $tableau_object;
}

sub random_float_matrix {

    # code to produce a matrix of random floats (or naturals)
    my $rows    = shift;
    my $columns = shift;
    my $natural_numbers;
    $natural_numbers = 0 unless $natural_numbers = shift;
    my $matrix;
    for my $i ( 0 .. $rows - 1 ) {
        for my $j ( 0 .. $columns - 1 ) {
            $matrix->[$i]->[$j] =
              $natural_numbers == 0 ? rand : int( 10 * rand );
        }
    }

    return $matrix;
}

sub matrix_copy {

    # code to copy matrix
    my $matrix = shift;
    my $matrix_copy;

    for my $i ( 0 .. $rows - 1 ) {
        for my $j ( 0 .. $columns - 1 ) {
            $matrix_copy->[$i]->[$j] = $matrix->[$i]->[$j];
        }
    }

    return $matrix_copy;
}
