package Algorithm::Simplex::Types;
use Moose::Util::TypeConstraints;
use PDL::Lite;
use Math::Cephes::Fraction qw(:fract);
use Data::Dumper;

# Want to coerce Array input into piddle.
subtype 'Piddle' 
    => as 'PDL' 
    => where { $_->isa('PDL') } 
    => message { "This thingy $_ is not a Piddle!" };

coerce 'Piddle' 
    => from 'ArrayRef[ArrayRef[Num]]' 
    => via { PDL->pdl($_) };

subtype 'FractMatrix' 
    => as 'ArrayRef[ArrayRef[Math::Cephes::Fraction]]' 
    => where { 1 } 
    => message { "This thingy $_ is not a matrix of Fraction objects.  It is a " . Dumper $_ };

coerce 'FractMatrix' 
    => from 'ArrayRef[ArrayRef[Num]]' 
    => via { &fraction_maker($_) };

=head2 fraction_maker

Make each integer and rational entry a Math::Cephes::Fraction object.

=cut

sub fraction_maker {
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

1;

