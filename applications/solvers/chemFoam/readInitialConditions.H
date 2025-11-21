
    const label nSpecie = Y.size();
    const scalarList W(::W(thermo));

    scalarList Y0(nSpecie, 0.0);
    scalarList X0(nSpecie, 0.0);

    dictionary fractions(initialConditions.subDict("fractions"));
    if (fractionBasis == "mole")
    {
        forAll(Y, i)
        {
            const word& name = Y[i].name();
            if (fractions.found(name))
            {
                X0[i] = fractions.lookup<scalar>(name);
            }
        }

        scalar mw = 0.0;
        const scalar mTot = sum(X0);
        forAll(Y, i)
        {
            X0[i] /= mTot;
            mw += W[i]*X0[i];
        }

        forAll(Y, i)
        {
            Y0[i] = X0[i]*W[i]/mw;
        }
    }
    else  // mass fraction
    {
        forAll(Y, i)
        {
            const word& name = Y[i].name();
            if (fractions.found(name))
            {
                Y0[i] = fractions.lookup<scalar>(name);
            }
        }

        scalar invW = 0.0;
        const scalar mTot = sum(Y0);
        forAll(Y, i)
        {
            Y0[i] /= mTot;
            invW += Y0[i]/W[i];
        }
        const scalar mw = 1.0/invW;

        forAll(Y, i)
        {
            X0[i] = Y0[i]*mw/W[i];
        }
    }

    // const scalar h0 = ::h0(thermo, Y0, p[0], T0);

    // forAll(Y, i)
    // {
    //     Y[i] = Y0[i];
    // }

    // thermo.he() = dimensionedScalar(dimEnergy/dimMass, h0);
    // thermo.correct();
