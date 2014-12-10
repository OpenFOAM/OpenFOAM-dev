{
    // Calculation of the mean radius based on SMR rs. Coefficient
    // factorGamma depends on nExp.
    scalar factorGamma = 1.;
    scalar delta = rs/factorGamma;

    scalar minValue = min(d/2.0, 0.04*rs);
    scalar maxValue = rs*4.0;

    scalar range = maxValue - minValue;

    if (maxValue - minValue < SMALL)
    {
        minValue = d/20.0;
        maxValue = d;
    }

    scalar nExp = 3.5;
    FixedList<scalar, 100> rrd;

    scalar probFactorMin = exp(-pow(minValue/delta, nExp));
    scalar probFactorMax = exp(-pow(maxValue/delta, nExp));
    scalar probFactor = 1./(probFactorMin - probFactorMax);

    forAll(rrd, n)
    {
        scalar xx = minValue + range*n/100;
        rrd[n] = (probFactorMin - exp(-pow(xx/delta, nExp)))*probFactor;
    }

    label n = 0;
    bool found = false;
    scalar random = rndGen.sample01<scalar>();

    while (!found && (n<100))
    {
        if (rrd[n] > random)
        {
            found = true;
        }
        n++;

    }

    rNew = minValue + range*(n - 0.5)/100.0;
}
