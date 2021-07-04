//Greg Miller (gmiller@gregmiller.net) 2021
//https://www.celestialprogramming.com/
//Released as public domain

/*
Conversion of original EPP82B.f Fortran code to JavaScript.
The code below has been reorganized and simplified from
the original Fortran code to make it easier to follow.  

The conversion and restructuring causes different rounding
for floating point.  Years near 2000 produce exactly the
same results, years 2000-3000 in the future and past differ
by 1e-6 which is an error of 1 micrometer, still well
within the margin of error.
*/

'use strict';

/*
Input :
tjj    julian date TDB (real double precision).
Output :
r[3]   table of rectangular coordinates (real double precision).
    reference frame : mean dynamical ecliptic and inertial
    equinox of J2000 (JD 2451545.0).
    r[1] : X (kilometer).
    r[2] : Y (kilometer).
    r[3] : Z (kilometer).
*/
class ELP82B{
    static ELP82B(tjj) {

        let t =  (tjj - 2451545.0) / 36525.;

        let r = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
            this.mainproblem(elp82data[i], i + 1, t, r);
        }

        for (let i = 3; i < 9; i++) {
            this.tides(elp82data[i], t, i + 1, r);
        }

        for (let i = 9; i < 21; i++) {
            this.planetary(elp82data[i], t, i + 1, r);
        }

        for (let i = 21; i < 36; i++) {
            this.tides(elp82data[i], t, i + 1, r);
        }

        //Change of coordinates.
        this.convertToJ2000(r, t);

        return r;
    }

    static convertToJ2000(r, t) {
        const rad = 648000. / Math.PI;
        const a0 = 384747.9806448954;
        const ath = 384747.9806743165;

        r[0] = r[0] / rad + 3.810344430588308 + 8399.684731773914*t + -0.000028547283984772807*t*t + 3.201709550047375e-8*t*t*t + -1.5363745554361197e-10*t*t*t*t
        r[1] = r[1] / rad
        r[2] = r[2] * a0 / ath
        
        let x1 = r[2] * Math.cos(r[1])
        const x2 = x1 * Math.sin(r[0])
        x1 = x1 * Math.cos(r[0])
        const x3 = r[2] * Math.sin(r[1])
        
        //Precession matrix.
        let pw = (0.10180391e-4 + 0.47020439e-6 * t + -0.5417367e-9 * t*t + -0.2507948e-11 * t*t*t + 0.463486e-14 * t*t*t*t) * t
        let qw = (-0.113469002e-3 + 0.12372674e-6 * t + 0.1265417e-8 * t*t + -0.1371808e-11 * t*t*t + -0.320334e-14 * t*t*t*t) * t
        
        const ra = 2. * Math.sqrt(1 - pw * pw - qw * qw)
        const pwqw = 2. * pw * qw
        const pw2 = 1 - 2. * pw * pw
        const qw2 = 1 - 2. * qw * qw
        pw = pw * ra
        qw = qw * ra
        r[0] = pw2 * x1 + pwqw * x2 + pw * x3
        r[1] = pwqw * x1 + qw2 * x2 - qw * x3
        r[2] = -pw * x1 + qw * x2 + (pw2 + qw2 - 1) * x3
    }

    //Main problem.
    static mainproblem(coeficients, ific, t, r) {

        const d = 5.198466741027443 + (7771.377146811758 + (-0.00002844935162118868 + (3.1973462269173895e-8 + -1.5436467606527627e-10*t)*t)*t)*t;
        const lp = -0.04312518020812495 + (628.301955168488 + (-0.000002680534842854624 + (7.126761112310179e-10 + 7.272205216643039e-13*t)*t)*t)*t;
        const l =  2.3555558982657994 + (8328.691426955555 + (0.00015702775761561094 + (2.504111144298864e-7 + -1.1863390776750345e-9*t)*t)*t)*t;
        const f = 1.627905233371468 + (8433.466158130539 + (-0.0000593921000043237 + (-4.949947684128362e-9 + 2.021673050226763e-11*t)*t)*t)*t;

        const rad = 648000. / Math.PI;
        const w12 = 1732559343.73604 / rad;

        const am = 0.074801329518;
        const dtasm = 2. * 0.002571881335 / (3. * am);
        
        const delnu = +0.55604 / rad / w12;
        const dele = +0.01789 / rad;
        const delg = -0.08066 / rad;
        const delnp = -0.06424 / rad / w12;
        const delep = -0.12879 / rad;

        const iv = ific - 1
        for (let i = 0; i < coeficients.length; i++) {
            const ilu = coeficients[i][0];
            const coef = coeficients[i][1];
            
            const tgv = coef[1] + dtasm * coef[5]
            let c0=coef[0]
            if (ific == 3) c0 = c0 - 2. * c0 * delnu / 3.
            const x = c0 + tgv * (delnp - am * delnu) + coef[2] * delg + coef[3] * dele + coef[4] * delep

            let y = ilu[0] * d + ilu[1] * lp + ilu[2] * l + ilu[3] * f;

            if (iv == 2) y = y + Math.PI / 2;
            y = y % (Math.PI * 2)
            r[iv] += x * Math.sin(y);
        }
    }

    //Planetary perturbations.
    static planetary(coeficients, t, ific, r) {

        const iv = ((ific - 1) % 3)

        const d = 5.198466741027443 + 7771.377146811758*t;
        const lp = -0.04312518020812495 + 628.301955168488*t;
        const l =  2.3555558982657994 + 8328.691426955555*t;
        const f = 1.627905233371468 + 8433.466158130539*t;

        const p1=(4.4026088424029615 + 2608.7903141574106 * t);
        const p2=(3.1761466969075944 + 1021.3285546211089 * t);
        const p3=(1.753470343150658 + 628.3075849621554 * t);
        const p4=(6.203480913399945 + 334.06124314922965 * t);
        const p5=(0.5995464973886735 + 52.96909650947205 * t);
        const p6=(0.8740167565184808 + 21.329909543800007 * t);
        const p7=(5.481293871604991 + 7.4781598567143535 * t);
        const p8=(5.311886286783467 + 3.813303563758456 * t);

        for (let i = 0; i < coeficients.length; i++) {
            const ipla = coeficients[i][0];
            const pha = coeficients[i][1];
            let x = coeficients[i][2];

            if (ific >= 13 && ific <= 15) x = x * t
            if (ific >= 19 && ific <= 21) x = x * t
            let y = pha * Math.PI / 180.0 + ipla[0] * p1 + ipla[1] * p2 + ipla[2] * p3 + ipla[3] * p4 + ipla[4] * p5 + ipla[5] * p6 + ipla[6] * p7;

            if (ific < 16) {
                y+=ipla[7] * p8 + ipla[8] * d + ipla[9] * l + ipla[10] * f;
            } else {
                y+=ipla[7] * d + ipla[8] * lp + ipla[9] * l + ipla[10] * f;
            }
            y = y % (Math.PI * 2)
            r[iv] += x * Math.sin(y)
        }
    }

    //Figures - Tides - Relativity - Solar eccentricity.
    static tides(coeficients, t, ific, r) {
        const iv = ((ific - 1) % 3)

        const d = 5.198466741027443 + 7771.377146811758*t;
        const lp = -0.04312518020812495 + 628.301955168488*t;
        const l =  2.3555558982657994 + 8328.691426955555*t;
        const f = 1.627905233371468 + 8433.466158130539*t;

        for (let j = 0; j < coeficients.length; j++) {
            const iz = coeficients[j][0];
            const ilu = coeficients[j][1];
            const pha = coeficients[j][2];
            let x = coeficients[j][3];

            if (ific >= 7 && ific <= 9) x = x * t
            if (ific >= 25 && ific <= 27) x = x * t
            if (ific >= 34 && ific <= 36) x = x * t*t

            let y = pha * Math.PI / 180.0 + ilu[0]*d + ilu[1]*lp + ilu[2]*l + ilu[3]*f + iz * (3.810344430588308 + 8399.709113522267*t);

            y = y % (Math.PI * 2);
            r[iv] += x * Math.sin(y)
        }
    }
}