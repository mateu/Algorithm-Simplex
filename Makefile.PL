use inc::Module::Install;

name 'Algorithm-Simplex';

all_from 'lib/Algorithm/Simplex.pm';

abstract 'Simplex Algorithm Implementation using Tucker Tableau';
author 'Mateu X Hunter';
license 'perl';
version_from 'lib/Algorithm/Simplex.pm';

requires(
    'Moose'                  => 0,
    'namespace::autoclean'   => 0,
    'Math::Cephes::Fraction' => 0,
    'PDL'                    => 0,
);

install_script glob('script/*.pl');
WriteAll;
