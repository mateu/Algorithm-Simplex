use inc::Module::Install;

name 'Algorithm-Simplex';

all_from 'lib/Algorithm/Simplex.pm';

abstract 'Simplex Algorithm Implementation using Tucker Tableau';
author 'Mateu X Hunter';
license 'perl';
version_from 'lib/Algorithm/Simplex.pm';


requires(
'Class::MethodMaker'            => 0,
'Math::Cephes::Fraction'        => 0,
'PDL'                           => 0,
'Readonly'                      => 0,
);


install_script glob('script/*.pl');
WriteAll;
