package Algorithm::Simplex::Types;
use Moose::Util::TypeConstraints;
use PDL::Lite;
use Math::Cephes::Fraction qw(:fract);
use Data::Dumper;

=head1 Name

Algorithm::Simplex::Types - Types into which we coerce matrix input
for PDL and Rational models

=cut

subtype 'Piddle' 
    => as 'PDL' 
    => where { $_->isa('PDL') }
    => message { "This thingy $_ is not a Piddle!" };
coerce 'Piddle' 
    => from 'ArrayRef[ArrayRef[Num]]' 
    => via { PDL->pdl($_) };
    
subtype 'PiddleDisplay'
    => as 'ArrayRef[ArrayRef[Str]]',
    => where { 1 }
    => message { "This thingys $_ is not a PiddleDisplay!"};
coerce 'PiddleDisplay'
    => from 'Piddle'
    => via { &display_piddle($_) };      

subtype 'FractionMatrix' 
    => as 'ArrayRef[ArrayRef[Math::Cephes::Fraction]]' 
    => where { 1 } 
    => message { "This thingy $_ is not a matrix of Fraction objects.  It is a " . Dumper $_ };
coerce 'FractionMatrix' 
    => from 'ArrayRef[ArrayRef[Num]]' 
    => via { &make_fractions($_) };
    
subtype 'FractionDisplay'
    => as 'ArrayRef[ArrayRef[Str]]'
    => where { 1 }
    => message { "This thingy $_ is not a FractDisplay!" };
coerce 'FractionDisplay'
    => from 'FractionMatrix'    
    => via { &display_fractions($_) };

=head1 Methods

=head2 make_fractions

Make each rational entry a Math::Cephes::Fraction object.

=cut

sub make_fractions {
    my $tableau = shift;
    
    for my $i ( 0 ..  scalar @{ $tableau } - 1 ) {
        for my $j ( 0 .. scalar @{ $tableau->[0] } - 1 ) {

            # Check for existing rationals indicated with "/"
            if ( $tableau->[$i]->[$j] =~ m{(\-?\d+)\/(\-?\d+)} ) {
                $tableau->[$i]->[$j] = fract( $1, $2 );
            }
            else {
                $tableau->[$i]->[$j] =
                  fract( $tableau->[$i]->[$j], 1 );
            }
        }
    }
    return $tableau;
}

=head display_fractions

Convert each fraction object entry into a string.

=cut

sub display_fractions {
    my $fraction_tableau = shift;

    my $display_tableau;
    for my $i ( 0 .. scalar @{ $fraction_tableau } - 1 ) {
        for my $j ( 0 .. scalar @{ $fraction_tableau->[0] } - 1  ) {
            $display_tableau->[$i]->[$j] = $fraction_tableau->[$i]->[$j]->as_string;
        }
    }
    return $display_tableau;

}

=head2 display_piddle 

Convert a PDL into an ArrayRef[ArrayRef[Num]]

=cut

sub display_piddle {
    my $piddle_tableau = shift;

    my @display_tableau;
    my ($number_of_columns, $number_of_rows) = ($piddle_tableau->dims);
    my $number_of_zero_based_rows = $number_of_rows - 1;
    my $number_of_zero_based_columns = $number_of_columns - 1;
    for my $i ( 0 .. $number_of_zero_based_rows ) {
        my $row = $piddle_tableau->slice("0:$number_of_zero_based_columns,($i)");
        my @row   = $row->list;
        push @display_tableau, \@row;
    }
    
    return \@display_tableau;
}

1;

