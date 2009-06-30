package Algorithm::Simplex::Rational;
use strict;
use warnings;
use base 'Algorithm::Simplex';
use Math::Cephes::Fraction qw(:fract);

my $one     = fract( 1, 1 );
my $neg_one = fract( 1, -1 );

=head1 Name

Algorithm::Simplex::Rational - Rational version of the Simplex Algorithm

=head1 Methods

=head2 pivot

Do the algebra of a Tucker/Bland pivot.  i.e. Traverse from one node to and 
adjacent node along the Simplex of feasible solutions.

=cut

sub pivot {

    my $self                = shift;
    my $pivot_row_number    = shift;
    my $pivot_column_number = shift;

    # Do tucker algebra on pivot row
    my $scale =
      $one->rdiv(
        $self->{_tableau}->[$pivot_row_number]->[$pivot_column_number] );
    for my $j ( 0 .. $self->{_number_of_columns} ) {
        $self->{_tableau}->[$pivot_row_number]->[$j] =
          $self->{_tableau}->[$pivot_row_number]->[$j]->rmul($scale);
    }
    $self->{_tableau}->[$pivot_row_number]->[$pivot_column_number] = $scale;

    # Do tucker algebra elsewhere
    for my $i ( 0 .. $self->{_number_of_rows} ) {
        if ( $i != $pivot_row_number ) {

            my $neg_a_ic =
              $self->{_tableau}->[$i]->[$pivot_column_number]->rmul($neg_one);
            for my $j ( 0 .. $self->{_number_of_columns} ) {
                $self->{_tableau}->[$i]->[$j] =
                  $self->{_tableau}->[$i]->[$j]->radd(
                    $neg_a_ic->rmul(
                        $self->{_tableau}->[$pivot_row_number]->[$j]
                    )
                  );
            }
            $self->{_tableau}->[$i]->[$pivot_column_number] =
              $neg_a_ic->rmul($scale);
        }
    }
}

sub determine_simplex_pivot_columns {
    my $self = shift;

    my @simplex_pivot_column_numbers;
    for my $col_num ( 0 .. $self->{_number_of_columns} - 1 ) {
        my $bottom_row_fraction =
          $self->{_tableau}->[ $self->{_number_of_rows} ]->[$col_num];
        my $bottom_row_numeric =
          $bottom_row_fraction->{n} / $bottom_row_fraction->{d};
        if ( $bottom_row_numeric > 0 ) {
            push( @simplex_pivot_column_numbers, $col_num );
        }
    }
    return (@simplex_pivot_column_numbers);
}

sub determine_positive_ratios {
    my $self                = shift;
    my $pivot_column_number = shift;

# Build Ratios and Choose row(s) that yields min for the bland simplex column as a candidate pivot point.
# To be a Simplex pivot we must not consider negative entries
    my %pivot_for;
    my @positive_ratios;
    my @positive_ratio_row_numbers;

    #print "Column: $possible_pivot_column\n";
    for my $row_num ( 0 .. $self->{_number_of_rows} - 1 ) {
        my $bottom_row_fraction =
          $self->{_tableau}->[$row_num]->[$pivot_column_number];
        my $bottom_row_numeric =
          $bottom_row_fraction->{n} / $bottom_row_fraction->{d};

        if ( $bottom_row_numeric > 0 ) {
            push(
                @positive_ratios,
                (
                    $self->{_tableau}->[$row_num]
                      ->[ $self->{_number_of_columns} ]->{n} *
                      $self->{_tableau}->[$row_num]->[$pivot_column_number]
                      ->{d}
                  ) / (
                    $self->{_tableau}->[$row_num]->[$pivot_column_number]
                      ->{n} *
                      $self->{_tableau}->[$row_num]
                      ->[ $self->{_number_of_columns} ]->{d}
                  )
            );

            # Track the rows that give ratios
            push @positive_ratio_row_numbers, $row_num;
        }
    }
    return ( \@positive_ratios, \@positive_ratio_row_numbers );
}
sub tableau_is_optimal {
    my $self = shift;

    # check basement row for having non-positive entries which
    # would => optimal when in phase 2.
    my $optimal_flag = 1;

    # if a positve entry exists in the basement row we don't have optimality
    for my $j ( 0 .. $self->{_number_of_columns} - 1 ) {
        my $basement_row_fraction = $self->{_tableau}->[ $self->{_number_of_rows} ]->[$j];
        my $basement_row_numeric = $basement_row_fraction->{n} / $basement_row_fraction->{d};
        if ( $basement_row_numeric > 0 ) {
            $optimal_flag = 0;
            last;
        }
    }
    return $optimal_flag;
}

sub convert_natural_number_tableau_to_fractional_object_tableau {
    my $self = shift;

# Make each integer and rational entry a fractional object for rational arthimetic
    for my $i ( 0 .. $self->{_number_of_rows} ) {
        for my $j ( 0 .. $self->{_number_of_columns} ) {

            # Check for existing rationals indicated with "/"
            if ( $self->{_tableau}->[$i]->[$j] =~ m{(\-?\d+)\/(\-?\d+)} ) {
                $self->{_tableau}->[$i]->[$j] = fract( $1, $2 );
            }
            else {
                $self->{_tableau}->[$i]->[$j] =
                  fract( $self->{_tableau}->[$i]->[$j], 1 );
            }
        }
    }
}

1;
