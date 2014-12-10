    dtChem = chemistry.solve(runTime.deltaT().value());
    scalar Sh = chemistry.Sh()()[0]/rho[0];
    integratedHeat += Sh*runTime.deltaT().value();
