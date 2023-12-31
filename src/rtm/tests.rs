use super::{liquid_cloud::*, oxygen::*, water_vapor::*, *};

use approx::assert_relative_eq;
use num_complex::Complex32;

/// Check some values for the Buck equation. These are compared against the
/// Fortran version.
#[test]
fn buck_vapor() {
    let input_temperatures = [273.15, 290., 300., 310.];
    let expected_outputs = [6.1121, 19.1925564, 35.3524513, 62.2872276];

    for (&temperature, &expected_output) in input_temperatures.iter().zip(&expected_outputs) {
        assert_relative_eq!(buck_vap(temperature), expected_output);
    }
}

/// Check some values for the water vapor absorption coefficient. These values
/// are from the Fortran version.
#[test]
fn water_vapor() {
    let inputs_and_outputs = [
        [
            1000.00000,
            273.149994,
            1.00000000,
            2.00000000,
            2.59790759E-05,
        ],
        [
            1000.00000,
            273.149994,
            1.00000000,
            10.0000000,
            7.06406950E-04,
        ],
        [
            1000.00000,
            273.149994,
            1.00000000,
            30.0000000,
            7.41568720E-03,
        ],
        [
            1000.00000,
            273.149994,
            0.100000001,
            2.00000000,
            2.57373767E-06,
        ],
        [
            1000.00000,
            273.149994,
            0.100000001,
            10.0000000,
            6.98733420E-05,
        ],
        [
            1000.00000,
            273.149994,
            0.100000001,
            30.0000000,
            7.33381894E-04,
        ],
        [
            1000.00000,
            273.149994,
            1.00000005E-03,
            2.00000000,
            2.57107917E-08,
        ],
        [
            1000.00000,
            273.149994,
            1.00000005E-03,
            10.0000000,
            6.97889163E-07,
        ],
        [
            1000.00000,
            273.149994,
            1.00000005E-03,
            30.0000000,
            7.32481112E-06,
        ],
        [
            1000.00000,
            290.000000,
            1.00000000,
            2.00000000,
            2.21492974E-05,
        ],
        [
            1000.00000,
            290.000000,
            1.00000000,
            10.0000000,
            6.05081732E-04,
        ],
        [
            1000.00000,
            290.000000,
            1.00000000,
            30.0000000,
            6.41573453E-03,
        ],
        [
            1000.00000,
            290.000000,
            0.100000001,
            2.00000000,
            2.19888739E-06,
        ],
        [
            1000.00000,
            290.000000,
            0.100000001,
            10.0000000,
            5.99988671E-05,
        ],
        [
            1000.00000,
            290.000000,
            0.100000001,
            30.0000000,
            6.36112352E-04,
        ],
        [
            1000.00000,
            290.000000,
            1.00000005E-03,
            2.00000000,
            2.19712266E-08,
        ],
        [
            1000.00000,
            290.000000,
            1.00000005E-03,
            10.0000000,
            5.99428290E-07,
        ],
        [
            1000.00000,
            290.000000,
            1.00000005E-03,
            30.0000000,
            6.35511469E-06,
        ],
        [
            1000.00000,
            300.000000,
            1.00000000,
            2.00000000,
            2.02244992E-05,
        ],
        [
            1000.00000,
            300.000000,
            1.00000000,
            10.0000000,
            5.53888734E-04,
        ],
        [
            1000.00000,
            300.000000,
            1.00000000,
            30.0000000,
            5.90435695E-03,
        ],
        [
            1000.00000,
            300.000000,
            0.100000001,
            2.00000000,
            2.00966224E-06,
        ],
        [
            1000.00000,
            300.000000,
            0.100000001,
            10.0000000,
            5.49828328E-05,
        ],
        [
            1000.00000,
            300.000000,
            0.100000001,
            30.0000000,
            5.86065173E-04,
        ],
        [
            1000.00000,
            300.000000,
            1.00000005E-03,
            2.00000000,
            2.00825561E-08,
        ],
        [
            1000.00000,
            300.000000,
            1.00000005E-03,
            10.0000000,
            5.49381696E-07,
        ],
        [
            1000.00000,
            300.000000,
            1.00000005E-03,
            30.0000000,
            5.85584257E-06,
        ],
        [
            500.000000,
            273.149994,
            1.00000000,
            2.00000000,
            1.31258521E-05,
        ],
        [
            500.000000,
            273.149994,
            1.00000000,
            10.0000000,
            3.61695827E-04,
        ],
        [
            500.000000,
            273.149994,
            1.00000000,
            30.0000000,
            3.94016085E-03,
        ],
        [
            500.000000,
            273.149994,
            0.100000001,
            2.00000000,
            1.28841123E-06,
        ],
        [
            500.000000,
            273.149994,
            0.100000001,
            10.0000000,
            3.53940814E-05,
        ],
        [
            500.000000,
            273.149994,
            0.100000001,
            30.0000000,
            3.85494815E-04,
        ],
        [
            500.000000,
            273.149994,
            1.00000005E-03,
            2.00000000,
            1.28575204E-08,
        ],
        [
            500.000000,
            273.149994,
            1.00000005E-03,
            10.0000000,
            3.53087728E-07,
        ],
        [
            500.000000,
            273.149994,
            1.00000005E-03,
            30.0000000,
            3.84557325E-06,
        ],
        [
            500.000000,
            290.000000,
            1.00000000,
            2.00000000,
            1.11654444E-05,
        ],
        [
            500.000000,
            290.000000,
            1.00000000,
            10.0000000,
            3.08859831E-04,
        ],
        [
            500.000000,
            290.000000,
            1.00000000,
            30.0000000,
            3.39363888E-03,
        ],
        [
            500.000000,
            290.000000,
            0.100000001,
            2.00000000,
            1.10049859E-06,
        ],
        [
            500.000000,
            290.000000,
            0.100000001,
            10.0000000,
            3.03698980E-05,
        ],
        [
            500.000000,
            290.000000,
            0.100000001,
            30.0000000,
            3.33620585E-04,
        ],
        [
            500.000000,
            290.000000,
            1.00000005E-03,
            2.00000000,
            1.09873346E-08,
        ],
        [
            500.000000,
            290.000000,
            1.00000005E-03,
            10.0000000,
            3.03131259E-07,
        ],
        [
            500.000000,
            290.000000,
            1.00000005E-03,
            30.0000000,
            3.32988748E-06,
        ],
        [
            500.000000,
            300.000000,
            1.00000000,
            2.00000000,
            1.01847809E-05,
        ],
        [
            500.000000,
            300.000000,
            1.00000000,
            10.0000000,
            2.82317982E-04,
        ],
        [
            500.000000,
            300.000000,
            1.00000000,
            30.0000000,
            3.11595527E-03,
        ],
        [
            500.000000,
            300.000000,
            0.100000001,
            2.00000000,
            1.00568730E-06,
        ],
        [
            500.000000,
            300.000000,
            0.100000001,
            10.0000000,
            2.78196840E-05,
        ],
        [
            500.000000,
            300.000000,
            0.100000001,
            30.0000000,
            3.06969858E-04,
        ],
        [
            500.000000,
            300.000000,
            1.00000005E-03,
            2.00000000,
            1.00428030E-08,
        ],
        [
            500.000000,
            300.000000,
            1.00000005E-03,
            10.0000000,
            2.77743482E-07,
        ],
        [
            500.000000,
            300.000000,
            1.00000005E-03,
            30.0000000,
            3.06460925E-06,
        ],
        [
            100.000000,
            273.149994,
            1.00000000,
            2.00000000,
            2.84019279E-06,
        ],
        [
            100.000000,
            273.149994,
            1.00000000,
            10.0000000,
            7.95400047E-05,
        ],
        [
            100.000000,
            273.149994,
            1.00000000,
            30.0000000,
            8.78573628E-04,
        ],
        [
            100.000000,
            273.149994,
            0.100000001,
            2.00000000,
            2.59843944E-07,
        ],
        [
            100.000000,
            273.149994,
            0.100000001,
            10.0000000,
            7.17569219E-06,
        ],
        [
            100.000000,
            273.149994,
            0.100000001,
            30.0000000,
            7.92029823E-05,
        ],
        [
            100.000000,
            273.149994,
            1.00000005E-03,
            2.00000000,
            2.57184651E-09,
        ],
        [
            100.000000,
            273.149994,
            1.00000005E-03,
            10.0000000,
            7.09007608E-08,
        ],
        [
            100.000000,
            273.149994,
            1.00000005E-03,
            30.0000000,
            7.82509801E-07,
        ],
        [
            100.000000,
            290.000000,
            1.00000000,
            2.00000000,
            2.37583617E-06,
        ],
        [
            100.000000,
            290.000000,
            1.00000000,
            10.0000000,
            6.66124397E-05,
        ],
        [
            100.000000,
            290.000000,
            1.00000000,
            30.0000000,
            7.41978583E-04,
        ],
        [
            100.000000,
            290.000000,
            0.100000001,
            2.00000000,
            2.21536595E-07,
        ],
        [
            100.000000,
            290.000000,
            0.100000001,
            10.0000000,
            6.14283135E-06,
        ],
        [
            100.000000,
            290.000000,
            0.100000001,
            30.0000000,
            6.83440303E-05,
        ],
        [
            100.000000,
            290.000000,
            1.00000005E-03,
            2.00000000,
            2.19771423E-09,
        ],
        [
            100.000000,
            290.000000,
            1.00000005E-03,
            10.0000000,
            6.08580493E-08,
        ],
        [
            100.000000,
            290.000000,
            1.00000005E-03,
            30.0000000,
            6.77000855E-07,
        ],
        [
            100.000000,
            300.000000,
            1.00000000,
            2.00000000,
            2.15075602E-06,
        ],
        [
            100.000000,
            300.000000,
            1.00000000,
            10.0000000,
            6.03525659E-05,
        ],
        [
            100.000000,
            300.000000,
            1.00000000,
            30.0000000,
            6.75204385E-04,
        ],
        [
            100.000000,
            300.000000,
            0.100000001,
            2.00000000,
            2.02283843E-07,
        ],
        [
            100.000000,
            300.000000,
            0.100000001,
            10.0000000,
            5.62106197E-06,
        ],
        [
            100.000000,
            300.000000,
            0.100000001,
            30.0000000,
            6.27955960E-05,
        ],
        [
            100.000000,
            300.000000,
            1.00000005E-03,
            2.00000000,
            2.00876737E-09,
        ],
        [
            100.000000,
            300.000000,
            1.00000005E-03,
            10.0000000,
            5.57549953E-08,
        ],
        [
            100.000000,
            300.000000,
            1.00000005E-03,
            30.0000000,
            6.22758421E-07,
        ],
    ];
    for [p, t, pv, freq, expected_output] in inputs_and_outputs {
        assert_relative_eq!(abh2o_rk_modified(p, t, pv, freq), expected_output);
    }
}

/// Check some values for the oxygen absorption coefficient. These values
/// are from the Fortran version.
#[test]
fn oxygen() {
    let inputs_and_outputs = [
        [
            1000.00000,
            273.149994,
            1.00000000,
            2.00000000,
            7.90619664E-03,
        ],
        [
            1000.00000,
            273.149994,
            1.00000000,
            10.0000000,
            9.75634996E-03,
        ],
        [
            1000.00000,
            273.149994,
            1.00000000,
            30.0000000,
            2.46610157E-02,
        ],
        [
            1000.00000,
            273.149994,
            0.100000001,
            2.00000000,
            7.91274104E-03,
        ],
        [
            1000.00000,
            273.149994,
            0.100000001,
            10.0000000,
            9.76428948E-03,
        ],
        [
            1000.00000,
            273.149994,
            0.100000001,
            30.0000000,
            2.46806331E-02,
        ],
        [
            1000.00000,
            273.149994,
            1.00000005E-03,
            2.00000000,
            7.91346002E-03,
        ],
        [
            1000.00000,
            273.149994,
            1.00000005E-03,
            10.0000000,
            9.76515841E-03,
        ],
        [
            1000.00000,
            273.149994,
            1.00000005E-03,
            30.0000000,
            2.46827919E-02,
        ],
        [
            1000.00000,
            290.000000,
            1.00000000,
            2.00000000,
            6.51325937E-03,
        ],
        [
            1000.00000,
            290.000000,
            1.00000000,
            10.0000000,
            7.95020070E-03,
        ],
        [
            1000.00000,
            290.000000,
            1.00000000,
            30.0000000,
            2.04622969E-02,
        ],
        [
            1000.00000,
            290.000000,
            0.100000001,
            2.00000000,
            6.51863590E-03,
        ],
        [
            1000.00000,
            290.000000,
            0.100000001,
            10.0000000,
            7.95667619E-03,
        ],
        [
            1000.00000,
            290.000000,
            0.100000001,
            30.0000000,
            2.04787478E-02,
        ],
        [
            1000.00000,
            290.000000,
            1.00000005E-03,
            2.00000000,
            6.51922543E-03,
        ],
        [
            1000.00000,
            290.000000,
            1.00000005E-03,
            10.0000000,
            7.95738958E-03,
        ],
        [
            1000.00000,
            290.000000,
            1.00000005E-03,
            30.0000000,
            2.04805583E-02,
        ],
        [
            1000.00000,
            300.000000,
            1.00000000,
            2.00000000,
            5.82994567E-03,
        ],
        [
            1000.00000,
            300.000000,
            1.00000000,
            10.0000000,
            7.07963528E-03,
        ],
        [
            1000.00000,
            300.000000,
            1.00000000,
            30.0000000,
            1.84053257E-02,
        ],
        [
            1000.00000,
            300.000000,
            0.100000001,
            2.00000000,
            5.83474990E-03,
        ],
        [
            1000.00000,
            300.000000,
            0.100000001,
            10.0000000,
            7.08540762E-03,
        ],
        [
            1000.00000,
            300.000000,
            0.100000001,
            30.0000000,
            1.84202138E-02,
        ],
        [
            1000.00000,
            300.000000,
            1.00000005E-03,
            2.00000000,
            5.83527796E-03,
        ],
        [
            1000.00000,
            300.000000,
            1.00000005E-03,
            10.0000000,
            7.08604185E-03,
        ],
        [
            1000.00000,
            300.000000,
            1.00000005E-03,
            30.0000000,
            1.84218492E-02,
        ],
        [
            500.000000,
            273.149994,
            1.00000000,
            2.00000000,
            2.12392258E-03,
        ],
        [
            500.000000,
            273.149994,
            1.00000000,
            10.0000000,
            2.44359416E-03,
        ],
        [
            500.000000,
            273.149994,
            1.00000000,
            30.0000000,
            6.15879428E-03,
        ],
        [
            500.000000,
            273.149994,
            0.100000001,
            2.00000000,
            2.12739035E-03,
        ],
        [
            500.000000,
            273.149994,
            0.100000001,
            10.0000000,
            2.44757137E-03,
        ],
        [
            500.000000,
            273.149994,
            0.100000001,
            30.0000000,
            6.16860390E-03,
        ],
        [
            500.000000,
            273.149994,
            1.00000005E-03,
            2.00000000,
            2.12777173E-03,
        ],
        [
            500.000000,
            273.149994,
            1.00000005E-03,
            10.0000000,
            2.44800886E-03,
        ],
        [
            500.000000,
            273.149994,
            1.00000005E-03,
            30.0000000,
            6.16968237E-03,
        ],
        [
            500.000000,
            290.000000,
            1.00000000,
            2.00000000,
            1.72994297E-03,
        ],
        [
            500.000000,
            290.000000,
            1.00000000,
            10.0000000,
            1.99030270E-03,
        ],
        [
            500.000000,
            290.000000,
            1.00000000,
            30.0000000,
            5.11018187E-03,
        ],
        [
            500.000000,
            290.000000,
            0.100000001,
            2.00000000,
            1.73276500E-03,
        ],
        [
            500.000000,
            290.000000,
            0.100000001,
            10.0000000,
            1.99354719E-03,
        ],
        [
            500.000000,
            290.000000,
            0.100000001,
            30.0000000,
            5.11840824E-03,
        ],
        [
            500.000000,
            290.000000,
            1.00000005E-03,
            2.00000000,
            1.73307571E-03,
        ],
        [
            500.000000,
            290.000000,
            1.00000005E-03,
            10.0000000,
            1.99390342E-03,
        ],
        [
            500.000000,
            290.000000,
            1.00000005E-03,
            30.0000000,
            5.11931209E-03,
        ],
        [
            500.000000,
            300.000000,
            1.00000000,
            2.00000000,
            1.53969473E-03,
        ],
        [
            500.000000,
            300.000000,
            1.00000000,
            10.0000000,
            1.77195808E-03,
        ],
        [
            500.000000,
            300.000000,
            1.00000000,
            30.0000000,
            4.59646620E-03,
        ],
        [
            500.000000,
            300.000000,
            0.100000001,
            2.00000000,
            1.54220534E-03,
        ],
        [
            500.000000,
            300.000000,
            0.100000001,
            10.0000000,
            1.77484925E-03,
        ],
        [
            500.000000,
            300.000000,
            0.100000001,
            30.0000000,
            4.60390979E-03,
        ],
        [
            500.000000,
            300.000000,
            1.00000005E-03,
            2.00000000,
            1.54248148E-03,
        ],
        [
            500.000000,
            300.000000,
            1.00000005E-03,
            10.0000000,
            1.77516695E-03,
        ],
        [
            500.000000,
            300.000000,
            1.00000005E-03,
            30.0000000,
            4.60472750E-03,
        ],
        [
            100.000000,
            273.149994,
            1.00000000,
            2.00000000,
            8.64337417E-05,
        ],
        [
            100.000000,
            273.149994,
            1.00000000,
            10.0000000,
            9.71217014E-05,
        ],
        [
            100.000000,
            273.149994,
            1.00000000,
            30.0000000,
            2.44592986E-04,
        ],
        [
            100.000000,
            273.149994,
            0.100000001,
            2.00000000,
            8.71413431E-05,
        ],
        [
            100.000000,
            273.149994,
            0.100000001,
            10.0000000,
            9.79184406E-05,
        ],
        [
            100.000000,
            273.149994,
            0.100000001,
            30.0000000,
            2.46556592E-04,
        ],
        [
            100.000000,
            273.149994,
            1.00000005E-03,
            2.00000000,
            8.72190940E-05,
        ],
        [
            100.000000,
            273.149994,
            1.00000005E-03,
            10.0000000,
            9.80059995E-05,
        ],
        [
            100.000000,
            273.149994,
            1.00000005E-03,
            30.0000000,
            2.46772397E-04,
        ],
        [
            100.000000,
            290.000000,
            1.00000000,
            2.00000000,
            7.01206882E-05,
        ],
        [
            100.000000,
            290.000000,
            1.00000000,
            10.0000000,
            7.90928752E-05,
        ],
        [
            100.000000,
            290.000000,
            1.00000000,
            30.0000000,
            2.02932788E-04,
        ],
        [
            100.000000,
            290.000000,
            0.100000001,
            2.00000000,
            7.06947685E-05,
        ],
        [
            100.000000,
            290.000000,
            0.100000001,
            10.0000000,
            7.97426983E-05,
        ],
        [
            100.000000,
            290.000000,
            0.100000001,
            30.0000000,
            2.04579439E-04,
        ],
        [
            100.000000,
            290.000000,
            1.00000005E-03,
            2.00000000,
            7.07578365E-05,
        ],
        [
            100.000000,
            290.000000,
            1.00000005E-03,
            10.0000000,
            7.98140973E-05,
        ],
        [
            100.000000,
            290.000000,
            1.00000005E-03,
            30.0000000,
            2.04760363E-04,
        ],
        [
            100.000000,
            300.000000,
            1.00000000,
            2.00000000,
            6.22867592E-05,
        ],
        [
            100.000000,
            300.000000,
            1.00000000,
            10.0000000,
            7.04105114E-05,
        ],
        [
            100.000000,
            300.000000,
            1.00000000,
            30.0000000,
            1.82524440E-04,
        ],
        [
            100.000000,
            300.000000,
            0.100000001,
            2.00000000,
            6.27967020E-05,
        ],
        [
            100.000000,
            300.000000,
            0.100000001,
            10.0000000,
            7.09894957E-05,
        ],
        [
            100.000000,
            300.000000,
            0.100000001,
            30.0000000,
            1.84014571E-04,
        ],
        [
            100.000000,
            300.000000,
            1.00000005E-03,
            2.00000000,
            6.28527414E-05,
        ],
        [
            100.000000,
            300.000000,
            1.00000005E-03,
            10.0000000,
            7.10531240E-05,
        ],
        [
            100.000000,
            300.000000,
            1.00000005E-03,
            30.0000000,
            1.84178294E-04,
        ],
    ];

    for [p, t, pv, freq, expected_output] in inputs_and_outputs {
        assert_relative_eq!(fdabsoxy_1992_modified(p, t, pv, freq), expected_output);
    }
}

/// Check some values for the liquid cloud absorption coefficient. These
/// values are from the Fortran version.
#[test]
fn liquid_cloud() {
    let inputs_and_outputs = [
        [2.00000000, 273.149994, 1.00000000, 8.60669592E-04],
        [2.00000000, 273.149994, 0.100000001, 8.60669606E-05],
        [2.00000000, 273.149994, 1.00000005E-03, 8.60669616E-07],
        [2.00000000, 290.000000, 1.00000000, 5.34837891E-04],
        [2.00000000, 290.000000, 0.100000001, 5.34837891E-05],
        [2.00000000, 290.000000, 1.00000005E-03, 5.34837909E-07],
        [2.00000000, 300.000000, 1.00000000, 4.27581428E-04],
        [2.00000000, 300.000000, 0.100000001, 4.27581472E-05],
        [2.00000000, 300.000000, 1.00000005E-03, 4.27581455E-07],
        [10.0000000, 273.149994, 1.00000000, 2.13160031E-02],
        [10.0000000, 273.149994, 0.100000001, 2.13160040E-03],
        [10.0000000, 273.149994, 1.00000005E-03, 2.13160038E-05],
        [10.0000000, 290.000000, 1.00000000, 1.33182406E-02],
        [10.0000000, 290.000000, 0.100000001, 1.33182411E-03],
        [10.0000000, 290.000000, 1.00000005E-03, 1.33182421E-05],
        [10.0000000, 300.000000, 1.00000000, 1.06627224E-02],
        [10.0000000, 300.000000, 0.100000001, 1.06627226E-03],
        [10.0000000, 300.000000, 1.00000005E-03, 1.06627231E-05],
        [30.0000000, 273.149994, 1.00000000, 0.178156421],
        [30.0000000, 273.149994, 0.100000001, 1.78156421E-02],
        [30.0000000, 273.149994, 1.00000005E-03, 1.78156421E-04],
        [30.0000000, 290.000000, 1.00000000, 0.116090909],
        [30.0000000, 290.000000, 0.100000001, 1.16090896E-02],
        [30.0000000, 290.000000, 1.00000005E-03, 1.16090901E-04],
        [30.0000000, 300.000000, 1.00000000, 9.40238163E-02],
        [30.0000000, 300.000000, 0.100000001, 9.40238219E-03],
        [30.0000000, 300.000000, 1.00000005E-03, 9.40238242E-05],
    ];

    for [freq, t, rhol, expected_output] in inputs_and_outputs {
        assert_relative_eq!(fdcldabs(freq, t, rhol), expected_output);
    }
}

/// Check some values for the water dielectric value. These
/// values are from the Fortran version.
#[test]
fn water_dielectric() {
    let inputs_and_outputs = [
        (
            2.00000000,
            273.149994,
            0.00000000,
            Complex32::new(83.9792633, -17.5693436),
        ),
        (
            2.00000000,
            273.149994,
            15.0000000,
            Complex32::new(80.2050705, -28.2318840),
        ),
        (
            2.00000000,
            273.149994,
            30.0000000,
            Complex32::new(76.7605515, -37.6490784),
        ),
        (
            2.00000000,
            290.000000,
            0.00000000,
            Complex32::new(80.1225433, -9.69443417),
        ),
        (
            2.00000000,
            290.000000,
            15.0000000,
            Complex32::new(76.3306427, -27.7735748),
        ),
        (
            2.00000000,
            290.000000,
            30.0000000,
            Complex32::new(72.8786850, -43.5699654),
        ),
        (
            2.00000000,
            300.000000,
            0.00000000,
            Complex32::new(77.0278702, -7.13627577),
        ),
        (
            2.00000000,
            300.000000,
            15.0000000,
            Complex32::new(73.3600235, -29.7685146),
        ),
        (
            2.00000000,
            300.000000,
            30.0000000,
            Complex32::new(70.0197601, -49.4610291),
        ),
        (
            10.0000000,
            273.149994,
            0.00000000,
            Complex32::new(42.1181870, -40.8917809),
        ),
        (
            10.0000000,
            273.149994,
            15.0000000,
            Complex32::new(41.4164009, -41.5028992),
        ),
        (
            10.0000000,
            273.149994,
            30.0000000,
            Complex32::new(40.9500046, -41.8011780),
        ),
        (
            10.0000000,
            290.000000,
            0.00000000,
            Complex32::new(58.8581696, -34.6061859),
        ),
        (
            10.0000000,
            290.000000,
            15.0000000,
            Complex32::new(56.4688148, -36.5488663),
        ),
        (
            10.0000000,
            290.000000,
            30.0000000,
            Complex32::new(54.3929825, -38.0325890),
        ),
        (
            10.0000000,
            300.000000,
            0.00000000,
            Complex32::new(63.3488808, -28.8427391),
        ),
        (
            10.0000000,
            300.000000,
            15.0000000,
            Complex32::new(60.4814758, -31.9897060),
        ),
        (
            10.0000000,
            300.000000,
            30.0000000,
            Complex32::new(57.9399567, -34.5400085),
        ),
        (
            30.0000000,
            273.149994,
            0.00000000,
            Complex32::new(12.3748817, -22.6335545),
        ),
        (
            30.0000000,
            273.149994,
            15.0000000,
            Complex32::new(12.1546707, -23.0304813),
        ),
        (
            30.0000000,
            273.149994,
            30.0000000,
            Complex32::new(12.3004131, -23.3680649),
        ),
        (
            30.0000000,
            290.000000,
            0.00000000,
            Complex32::new(21.5122662, -30.7900391),
        ),
        (
            30.0000000,
            290.000000,
            15.0000000,
            Complex32::new(20.7252617, -30.8957348),
        ),
        (
            30.0000000,
            290.000000,
            30.0000000,
            Complex32::new(20.3074341, -30.8547306),
        ),
        (
            30.0000000,
            300.000000,
            0.00000000,
            Complex32::new(27.9100685, -33.4010315),
        ),
        (
            30.0000000,
            300.000000,
            15.0000000,
            Complex32::new(26.6922417, -33.5059967),
        ),
        (
            30.0000000,
            300.000000,
            30.0000000,
            Complex32::new(25.8418388, -33.4098701),
        ),
    ];

    for (freq, temp, s, expected_output) in inputs_and_outputs {
        let actual_output = meissner(freq, temp, s);
        // I'm not sure why, but the real part just needs a little higher
        // epsilon than the default, which is `f32::EPSILON`.
        assert_relative_eq!(actual_output.re, expected_output.re, epsilon = 5e-6);
        assert_relative_eq!(actual_output.im, expected_output.im);
    }
}
