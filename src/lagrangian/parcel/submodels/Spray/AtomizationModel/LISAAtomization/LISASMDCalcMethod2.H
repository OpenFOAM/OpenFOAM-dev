{
    // calculate the new diameter with the standard 1D Rosin Rammler
    // distribution. Calculation of the mean radius based on SMR rs.
    // Coefficient factorGamma depends on nExp. Note that Reitz either used
    // (Schmidt et al., 1999-01-0496) or skipped (Senecal et al.) this factor!
    // scalar factorGamma = 0.75*sqrt(mathematicalConstant::pi);   // nExp=2
    scalar factorGamma = 1.;
    scalar delta = dD/factorGamma;

    // dD is the SMD, and the delta is calculated using gamma
    // function. Here we assume nExp = 2
    scalar minValue = dD/10.0;
    scalar maxValue = dD;

    // The pdf value for 4.0*delta is already very small.
    // scalar maxValue = delta*4.0;

    if (maxValue - minValue < small)
    {
        minValue = maxValue/20.0;
    }

    scalar range = maxValue - minValue;

    scalar nExp = 3;
    FixedList<scalar, 500> rrd;
    scalar probFactorMin = exp(-pow(minValue/delta, nExp));
    scalar probFactorMax = exp(-pow(maxValue/delta, nExp));
    scalar probFactor = 1.0/(probFactorMin - probFactorMax);

    forAll(rrd, n)
    {
        scalar xx = minValue + range*n/500;
        rrd[n] = (probFactorMin - exp(-pow(xx/delta, nExp)))*probFactor;
    }


    bool success = false;

    scalar y = rndGen.sample01<scalar>();
    label k = 0;

    while(!success && (k<500))
    {
        if (rrd[k] > y)
        {
            success = true;
        }
        k++;
    }

    x = minValue + range*(k - 0.5)/500.0;
}
