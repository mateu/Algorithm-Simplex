use inc::Module::Install;

name 'Algorithm-Simplex';

all_from 'lib/Algorithm/Simplex.pm';

abstract 'Simplex Algorithm Implementation using Tucker Tableau';
author 'Mateu X Hunter';
license 'perl';
version_from 'lib/Algorithm/Simplex.pm';

requires 'Math::Cephes::Fraction';
requires 'PDL';
requires 'Class::MethodMaker';
requires 'HTML::Table';
requires 'HTML::TreeBuilder';
requires 'Readonly';

WriteAll;
