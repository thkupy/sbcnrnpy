TITLE na_pxko.mod  A fast voltage activated sodium channel with reduced inactivation

COMMENT
Created by TK (08-07-2016) for the pxko project
ENDCOMMENT

NEURON {
    SUFFIX na_pxko
    USEION na READ ena WRITE ina
    RANGE gnabar, ina, haP1, haP2, hbP1, hbP2
    GLOBAL minf, hinf, mtau, htau
    THREADSAFE
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gnabar=.120 (mho/cm2) <0,1e9>
    haP1=90.0
    haP2=90.0
    hbP1=30.0
    hbP2=15.0
}

STATE {
    m h
}

ASSIGNED {
    v (mV)
    celsius (degC) : 6.3
    ena (mV)
    ina (mA/cm2)
    minf hinf
    mtau (ms)
    htau (ms)
}

INITIAL {
    rates(v)
    m = minf
    h = hinf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gnabar*m*m*m*h*(v - ena)
}

DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
}

FUNCTION malpha(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3^((celsius - 22(degC))/10(degC))

    malpha = 0.36 * Tf * expM1( -(v+49), 3 )
    :malpha = 0.36 * Tf * (v + 49) / (1 - exp(-(v+49)/3))
}

FUNCTION mbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3^((celsius - 22(degC))/10(degC))

    mbeta = 0.4 * Tf * expM1( (v+58), 20 )
    :mbeta = -0.4 * Tf * (v + 58) / (1 - exp((v+58)/20))
}

FUNCTION halpha(v(mV)) (/ms) {
    LOCAL Tf, Tf10
    Tf = 3^((celsius - 22(degC))/10(degC))
    Tf10 = 10^((celsius - 22(degC))/10(degC))

    halpha = 2.4 * Tf / (1 + exp((v+haP1)/3.))  +  0.8 * Tf10 / (1 + exp(v+haP2))
}

FUNCTION hbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3^((celsius - 22(degC))/10(degC))

    hbeta = 3.6 * Tf / (1 + exp(-(v+hbP1)/hbP2))
}

FUNCTION expM1(x,y) {
    if (fabs(x/y) < 1e-6) {
	expM1 = y*(1 - x/y/2)
    } else {
	expM1 = x/(exp(x/y) - 1)
    }
}

PROCEDURE rates(v(mV)) {
    TABLE minf, hinf, mtau, htau DEPEND celsius FROM -100 TO 100 WITH 200

    mtau = 1/(malpha(v) + mbeta(v))
    minf = malpha(v)/(malpha(v) + mbeta(v))

    htau = 1/(halpha(v) + hbeta(v))
    hinf = halpha(v)/(halpha(v) + hbeta(v))
}
