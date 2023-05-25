if (constProp != "temperature")
{
    volScalarField& h = thermo.he();

    if (constProp == "volume")
    {
        h[0] = u0 + p[0]/rho[0] + integratedHeat;
    }
    else
    {
        h[0] = h0 + integratedHeat;
    }

    thermo.correct();
}
