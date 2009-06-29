package Algorithm::Simplex;
use strict;
use warnings;
use vars '$AUTOLOAD';    # Keep 'use strict' happy
use Carp;

=head1 Name

Algorithm::Simplex

=head1 Methods

=head2 new

Create a new tableau which is a ArrayRef[ArrayRef], i.e. two-dimensional array.

=cut

sub new {
    my $class = $_[0];
    bless { _tableau => $_[1], }, $class;
}

use Class::MethodMaker abstract => [
    qw(
      tableau_is_optimal
      pivot
      )
];

{
    my %_attrs = (
        _tableau           => 'read',
        _number_of_rows    => 'read',
        _number_of_columns => 'read',
    );

    sub _accessible {
        my ( $self, $attr, $mode ) = @_;
        $_attrs{$attr} =~ m{$mode};
    }
}

sub Algorithm::Simplex::AUTOLOAD {
    no strict "refs";
    my ( $self, $newval ) = @_;

    # Handle get_ methods
    if ( $AUTOLOAD =~ m{.*::get(_\w+)} && $self->_accessible( $1, 'read' ) ) {
        my $attr_name = $1;
        *{$AUTOLOAD} = sub { return $_[0]->{$attr_name} };
        return $self->{$attr_name};
    }

    # Otherwise a method has been called that doesn't exist
    croak "No such method: $AUTOLOAD";
}

sub Algorithm::Simplex::DESTROY {
    my $self = shift;
}

=head2 set_generic_variable_names_from_dimensions

Create variable names: x1, x2 ... , y1, y2, ... , u1, u2 ... , v1, v2 ...

Our variables are represented by:

    x, y, u, and v 
    
as found in Nering and Tuckers' book. 

x and y are for the primal LP while u and v belong to the dual LP.

=cut

sub set_generic_variable_names_from_dimensions {
    my $self = shift;
    my ( @x, @y, @v, @u );
    for my $i ( 0 .. $self->{_number_of_rows} - 1 ) {
        my $tmp_num = $i + 1;
        my $y       = 'y' . $tmp_num;
        $self->{_y_variables}->[$i]->{'generic'} = $y;
        my $v = 'v' . $tmp_num;
        $self->{_v_variables}->[$i]->{'generic'} = $v;
    }
    for my $j ( 0 .. $self->{_number_of_columns} - 1 ) {
        my $tmp_num = $j + 1;
        my $x       = 'x' . $tmp_num;
        $self->{_x_variables}->[$j]->{'generic'} = $x;
        my $u = 'u' . $tmp_num;
        $self->{_u_variables}->[$j]->{'generic'} = $u;
    }
}

=head2 get_bland_number_for

Given a column number (which represents a u variable) build the bland number 
from the generic variable name.

=cut

sub get_bland_number_for {
    my $self          = shift;
    my $variable_type = shift;
    my $variables     = '_' . $variable_type . '_variables';
    my $index         = shift;
    my $generic_name  = $self->{$variables}->[$index]->{'generic'};
    $generic_name =~ m{(.)(\d+)};
    my $var = $1;
    my $num = $2;
    my $start_num =
        $var eq 'x' ? 1
      : $var eq 'y' ? 2
      : $var eq 'v' ? 4
      : $var eq 'u' ? 3
      :               die "Variable name: $var does not equal x, y, v or u";
    my $bland_number = $start_num . $num;
    return $bland_number;
}

=head2 determine_bland_pivot_column_number

Find the pivot column using Bland ordering technique to prevent cycles.

=cut

sub determine_bland_pivot_column_number {
    my $self                         = shift;
    my @simplex_pivot_column_numbers = @_;

    my @bland_number_for_simplex_pivot_column;
    foreach my $col_number (@simplex_pivot_column_numbers) {
        push @bland_number_for_simplex_pivot_column,
          $self->get_bland_number_for( 'x', $col_number );
    }

# Pass blands number to routine that returns index of location where minimum bland occurs.
# Use this index to return the bland column column number from @positive_profit_column_numbers
    my @bland_column_number_index =
      $self->min_index( \@bland_number_for_simplex_pivot_column );
    my $bland_column_number_index = $bland_column_number_index[0];

    return $simplex_pivot_column_numbers[$bland_column_number_index];
}

=head2 determine_bland_pivot_row_number

Find the pivot row using Bland ordering technique to prevent cycles.

=cut

sub determine_bland_pivot_row_number {
    my $self = shift;
    my ( $positive_ratios, $positive_ratio_row_numbers ) = @_;

   # Now that we have the ratios and their respective rows we can find the min
   # and then select the lowest bland min if there are ties.
    my @min_indices = $self->min_index($positive_ratios);
    my @min_ratio_row_numbers =
      map { $positive_ratio_row_numbers->[$_] } @min_indices;
    my @bland_number_for_min_ratio_rows;
    foreach my $row_number (@min_ratio_row_numbers) {
        push @bland_number_for_min_ratio_rows,
          $self->get_bland_number_for( 'y', $row_number );
    }

# Pass blands number to routine that returns index of location where minimum bland occurs.
# Use this index to return the bland row number.
    my @bland_min_ratio_row_index =
      $self->min_index( \@bland_number_for_min_ratio_rows );
    my $bland_min_ratio_row_index = $bland_min_ratio_row_index[0];
    return $min_ratio_row_numbers[$bland_min_ratio_row_index];
}

=head2 min_index

Detemine the index of the element with minimal value.  
Used when finding bland pivots.

=cut

sub min_index {
    my $self = shift;
    my $l    = $_[0];
    my $n    = @{$l};
    return () unless $n;
    my $v_min = $l->[0];
    my @i_min = (0);

    for ( my $i = 1 ; $i < $n ; $i++ ) {
        if ( $l->[$i] < $v_min ) {
            $v_min = $l->[$i];
            @i_min = ($i);
        }
        elsif ( $l->[$i] == $v_min ) {
            push @i_min, $i;
        }
    }
    return @i_min;

}

=head2 exchange_pivot_variables

Exchange the variables when the a pivot is done.  The method pivot does the
algrebra while this method does the variable swapping.

=cut

sub exchange_pivot_variables {
    my $self                = shift;
    my $pivot_row_number    = shift;
    my $pivot_column_number = shift;

    # exchange variables based on $pivot_column_number and $pivot_row_number
    my $increasing_primal_variable =
      $self->{_x_variables}->[$pivot_column_number];
    my $zeroeing_primal_variable = $self->{_y_variables}->[$pivot_row_number];
    $self->{_x_variables}->[$pivot_column_number] = $zeroeing_primal_variable;
    $self->{_y_variables}->[$pivot_row_number] = $increasing_primal_variable;

    my $increasing_dual_variable = $self->{_v_variables}->[$pivot_row_number];
    my $zeroeing_dual_variable =
      $self->{_u_variables}->[$pivot_column_number];
    $self->{_v_variables}->[$pivot_row_number]    = $zeroeing_dual_variable;
    $self->{_u_variables}->[$pivot_column_number] = $increasing_dual_variable;
}

=head2 set_number_of_rows_and_columns

Given a tableau (matrix), determine its size.

=cut


sub set_number_of_rows_and_columns {
    my $self = shift;

    my @rows           = @{ $self->{_tableau} };
    my $number_of_rows = @rows;
    $number_of_rows -= 1;
    $self->{_number_of_rows} = $number_of_rows;

    my @columns           = @{ $self->{_tableau}->[0] };
    my $number_of_columns = @columns;
    $number_of_columns -= 1;
    $self->{_number_of_columns} = $number_of_columns;

}

=head2 get_row_and_column_numbers 

Get the dimensions of the tableau.

=cut

sub get_row_and_column_numbers {
    my $self = shift;
    return $self->{_number_of_rows}, $self->{_number_of_columns};
}

=head2 determine_bland_pivot_row_and_column_numbers

Higher level functions that uses other to return the (bland) pivot point.

=cut


sub determine_bland_pivot_row_and_column_numbers {
    my $self = shift;

    #return "you";
    my @simplex_pivot_columns = $self->determine_simplex_pivot_columns;
    my $tmp_out               = "pivot column: " . $simplex_pivot_columns[0];

    #return $tmp_out;
    my $pivot_column_number =
      $self->determine_bland_pivot_column_number(@simplex_pivot_columns);

    #return "pivot column number: $pivot_column_number";
    my ( $positive_ratios, $positive_ratio_row_numbers ) =
      $self->determine_positive_ratios($pivot_column_number);

    #return "postive_ratio: $positive_ratios->[0]";
    my $pivot_row_number =
      $self->determine_bland_pivot_row_number( $positive_ratios,
        $positive_ratio_row_numbers );

    #return "pivot row: $pivot_row_number";
    return ( $pivot_row_number, $pivot_column_number );
}

1;

__END__

=head1 Description

Base class for the Simplex model.  It defines some of the methods concretely,
and others such as:

=over 3

=item *

pivot

=item *

tableau_is_optimal

=item *

determine_positive_ratios

=item *

deteremine_simplex_pivot_columns

=back

are implemented in one of the three model types:

=over 3

=item *

Float

=item *

Rational

=item *

PDL

=back

=cut
