//Greg Miller (gmiller@gregmiller.net) 2021
//https://www.celestialprogramming.com/
//Released as public domain
//Conversion of original EPP82B Fortran code to JavaScript
'use strict';
function ELP82B(tjj) {
    /*
    //     Input :
    //     tjj    julian date TDB (real double precision).
    //
    //     Output :
    //     r[3]   table of rectangular coordinates (real double precision).
    //            reference frame : mean dynamical ecliptic and inertial
    //            equinox of J2000 (JD 2451545.0).
    //            r[1] : X (kilometer).
    //            r[2] : Y (kilometer).
    //            r[3] : Z (kilometer).
    */

    //     Initialisation.
    ///
    let r = [0, 0, 0];

    const deg = Math.PI / 180.;
    const rad = 648000. / Math.PI;
    const am = 0.074801329518;
    const alfa = 0.002571881335;
    const dtasm = 2. * alfa / (3. * am);
    //
    //     Lunar arguments.
    //
    let w = [[0], [0], [0], [0]];
    let eart = [];
    let peri = [];
    w[1][1] = (218 + 18 / 60.0 + 59.95571 / 3600.0) * deg;
    w[1][2] = 1732559343.73604 / rad;
    w[1][3] = -5.8883 / rad;
    w[1][4] = 0.6604e-2 / rad;
    w[1][5] = -0.3169e-4 / rad;
    w[2][1] = (83 + 21 / 60.0 + 11.67475 / 3600.0) * deg;
    w[2][2] = 14643420.2632 / rad;
    w[2][3] = -38.2776 / rad;
    w[2][4] = -0.45047e-1 / rad;
    w[2][5] = 0.21301e-3 / rad;
    w[3][1] = (125 + 2 / 60.0 + 40.39816 / 3600.0) * deg;
    w[3][2] = -6967919.3622 / rad;
    w[3][3] = 6.3622 / rad;
    w[3][4] = 0.7625e-2 / rad;
    w[3][5] = -0.3586e-4 / rad;
    eart[1] = (100 + 27 / 60.0 + 59.22059 / 3600.0) * deg;
    eart[2] = 129597742.2758 / rad;
    eart[3] = -0.0202 / rad;
    eart[4] = 0.9e-5 / rad;
    eart[5] = 0.15e-6 / rad;
    peri[1] = (102 + 56 / 60.0 + 14.42753 / 3600.0) * deg;
    peri[2] = 1161.2283 / rad;
    peri[3] = 0.5327 / rad;
    peri[4] = -0.138e-3 / rad;
    peri[5] = 0.;
    //
    //     Planetary arguments.
    //
    const preces = 5029.0966 / rad;

    //
    //     Delaunay's arguments.
    //
    let del = [[], [], [], [], [],];
    for (let i = 1; i <= 5; i++) {
        del[1][i] = w[1][i] - eart[i]; //d
        del[2][i] = eart[i] - peri[i]; //lp
        del[3][i] = w[1][i] - w[2][i]; //l
        del[4][i] = w[1][i] - w[3][i]; //f
    }
    del[1][1] = del[1][1] + Math.PI;

    let zeta = [];
    zeta[1] = w[1][1];
    zeta[2] = w[1][2] + preces;

    let t = [];
    t[1] = 1.
    t[2] = (tjj - 2451545.0) / 36525.;
    t[3] = t[2] * t[2];
    t[4] = t[3] * t[2];
    t[5] = t[4] * t[2];

    for (let i = 0; i < 3; i++) {
        mainproblem(data[i], t, del, dtasm, am, i + 1, r)
    }

    for (let i = 3; i < 9; i++) {
        tides(data[i], t, del, zeta, i + 1, r);
    }

    for (let i = 9; i < 21; i++) {
        planetary(data[i], t, del, i + 1, r);
    }

    for (let i = 21; i < 36; i++) {
        tides(data[i], t, del, zeta, i + 1, r);
    }

    //Change of coordinates.
    convertToJ2000(r, t, w);

    return r;
}

function convertToJ2000(r, t, w) {
    const rad = 648000. / Math.PI;
    const a0 = 384747.9806448954;
    const ath = 384747.9806743165;

    //     Precession matrix.
    const p1 = 0.10180391e-4;
    const p2 = 0.47020439e-6;
    const p3 = -0.5417367e-9;
    const p4 = -0.2507948e-11;
    const p5 = 0.463486e-14;
    const q1 = -0.113469002e-3;
    const q2 = 0.12372674e-6;
    const q3 = 0.1265417e-8;
    const q4 = -0.1371808e-11;
    const q5 = -0.320334e-14;

    r[0] = r[0] / rad + w[1][1] + w[1][2] * t[2] + w[1][3] * t[3] + w[1][4] * t[4] + w[1][5] * t[5]
    r[1] = r[1] / rad
    r[2] = r[2] * a0 / ath
    let x1 = r[2] * Math.cos(r[1])
    const x2 = x1 * Math.sin(r[0])
    x1 = x1 * Math.cos(r[0])
    const x3 = r[2] * Math.sin(r[1])
    let pw = (p1 + p2 * t[2] + p3 * t[3] + p4 * t[4] + p5 * t[5]) * t[2]
    let qw = (q1 + q2 * t[2] + q3 * t[3] + q4 * t[4] + q5 * t[5]) * t[2]
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
function mainproblem(coeficients, t, del, dtasm, am, ific, r) {
    const rad = 648000. / Math.PI;
    const w12 = 1732559343.73604 / rad;

    const delnu = +0.55604 / rad / w12;
    const dele = +0.01789 / rad;
    const delg = -0.08066 / rad;
    const delnp = -0.06424 / rad / w12;
    const delep = -0.12879 / rad;

    const iv = ((ific - 1) % 3) + 1
    for (let l = 0; l < coeficients.length; l++) {
        const ilu = coeficients[l][0];
        const coef = coeficients[l][1];
        const tgv = coef[1] + dtasm * coef[5]
        if (ific == 3) coef[0] = coef[0] - 2. * coef[0] * delnu / 3.
        const x = coef[0] + tgv * (delnp - am * delnu) + coef[2] * delg + coef[3] * dele + coef[4] * delep
        let y = 0.

        const T=t[1];
        const f  = (335779.55755 +
            (1739527263.0983 +
            (-12.2505 +
            (-0.001021 + 0.00000417 * T) * T) * T) * T) /rad;

        y = y + ilu[0] * (del[1][1] * t[1] + del[1][2] * t[2] + del[1][3] * t[3] + del[1][4] * t[4] + del[1][5] * t[5])
        y = y + ilu[1] * (del[2][1] * t[1] + del[2][2] * t[2] + del[2][3] * t[3] + del[2][4] * t[4] + del[2][5] * t[5])
        y = y + ilu[2] * (del[3][1] * t[1] + del[3][2] * t[2] + del[3][3] * t[3] + del[3][4] * t[4] + del[3][5] * t[5])
        //y = y + ilu[3] * f;
        y = y + ilu[3] * (del[4][1] * t[1] + del[4][2] * t[2] + del[4][3] * t[3] + del[4][4] * t[4] + del[4][5] * t[5])

        if(l==0){
            console.log(f,(del[4][5])*rad,-0.001021 + 0.00000417 * T);
        }
        
        if (iv == 3) y = y + Math.PI / 2;
        y = y % (Math.PI * 2);
        r[iv - 1] = r[iv - 1] + x * Math.sin(y);
    }
}

//Planetary perturbations.
function planetary(coeficients, t, del, ific, r) {
    const deg = Math.PI / 180.;
    const rad = 648000. / Math.PI;

    let eart = [];
    eart[1] = (100 + 27 / 60.0 + 59.22059 / 3600.0) * deg;
    eart[2] = 129597742.2758 / rad;

    let p = [[], [], [], [], [], [], [], [], []];

    p[1][1] = (252 + 15 / 60.0 + 3.25986 / 3600.0) * deg;
    p[2][1] = (181 + 58 / 60.0 + 47.28305 / 3600.0) * deg;
    p[3][1] = eart[1];
    p[4][1] = (355 + 25 / 60.0 + 59.78866 / 3600.0) * deg;
    p[5][1] = (34 + 21 / 60.0 + 5.34212 / 3600.0) * deg;
    p[6][1] = (50 + 4 / 60.0 + 38.89694 / 3600.0) * deg;
    p[7][1] = (314 + 3 / 60.0 + 18.01841 / 3600.0) * deg;
    p[8][1] = (304 + 20 / 60.0 + 55.19575 / 3600.0) * deg;
    p[1][2] = 538101628.68898 / rad;
    p[2][2] = 210664136.43355 / rad;
    p[3][2] = eart[2];
    p[4][2] = 68905077.59284 / rad;
    p[5][2] = 10925660.42861 / rad;
    p[6][2] = 4399609.65932 / rad;
    p[7][2] = 1542481.19393 / rad;
    p[8][2] = 786550.32074 / rad;

    const iv = ((ific - 1) % 3) + 1
    for (let l = 0; l < coeficients.length; l++) {
        const ipla = coeficients[l][0];
        const pha = coeficients[l][1];
        let x = coeficients[l][2];
        const per = coeficients[l][3];

        if (ific >= 13 && ific <= 15) x = x * t[2]
        if (ific >= 19 && ific <= 21) x = x * t[2]
        let y = pha * Math.PI / 180.0

        if (ific < 16) {
            for (let k = 1; k <= 2; k++) {
                y = y + (ipla[8] * del[1][k] + ipla[9] * del[3][k] + ipla[10] * del[4][k]) * t[k]
                for (let i = 1; i <= 8; i++) {
                    y = y + ipla[i - 1] * p[i][k] * t[k]
                }
            }
        } else {
            for (let k = 1; k <= 2; k++) {
                for (let i = 1; i <= 4; i++) {
                    y = y + ipla[i + 7 - 1] * del[i][k] * t[k]
                }
                for (let i = 1; i <= 7; i++) {
                    y = y + ipla[i - 1] * p[i][k] * t[k]
                }
            }
        }
        y = y % (Math.PI * 2)
        r[iv - 1] = r[iv - 1] + x * Math.sin(y)
    }
}

//Figures - Tides - Relativity - Solar eccentricity.
function tides(coeficients, t, del, zeta, ific, r) {
    const iv = ((ific - 1) % 3) + 1
    for (let l = 0; l < coeficients.length; l++) {
        const iz = coeficients[l][0];
        const ilu = coeficients[l][1];
        const pha = coeficients[l][2];
        let x = coeficients[l][3];
        const per = coeficients[l][4];

        if (ific >= 7 && ific <= 9) x = x * t[2]
        if (ific >= 25 && ific <= 27) x = x * t[2]
        if (ific >= 34 && ific <= 36) x = x * t[3]
        let y = pha * Math.PI / 180.0
        for (let k = 1; k <= 2; k++) {
            y = y + iz * zeta[k] * t[k]
            for (let i = 1; i <= 4; i++) {
                y = y + ilu[i - 1] * del[i][k] * t[k]
            }
        }
        y = y % (Math.PI * 2);
        r[iv - 1] = r[iv - 1] + x * Math.sin(y)
    }
}