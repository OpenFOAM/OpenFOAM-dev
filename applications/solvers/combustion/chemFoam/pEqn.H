{
    rho = thermo.rho();
    if (constProp == "volume")
    {
        scalar invW = 0.0;
        forAll(Y, i)
        {
            invW += Y[i][0]/W[i];
        }

        Rspecific[0] = 1000.0*constant::physicoChemical::R.value()*invW;

        p[0] = rho0*Rspecific[0]*thermo.T()[0];
        rho[0] = rho0;
    }
}
