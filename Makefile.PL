use inc::Module::Install;

name 'Algorithm-Simplex';

all_from 'lib/Algorithm/Simplex.pm';

abstract 'Simplex Algorithm Implementation using Tucker Tableaux';
author 'Mateu X Hunter';
license 'perl';
version_from 'lib/Algorithm/Simplex.pm';

requires(
    'Moo'                    => 0,
    'namespace::clean'       => 0,
    'Math::Cephes::Fraction' => 0,
    'PDL'                    => 0,
    'PDL::Lite'              => 0,
    'Math::BigRat'           => 0,
);

install_script glob('script/*.pl');
auto_install;
WriteAll;
