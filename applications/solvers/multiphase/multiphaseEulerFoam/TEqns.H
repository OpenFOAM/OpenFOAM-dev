{
    volScalarField kByCp1("kByCp1", alpha1*(k1/Cp1/rho1 + sqr(Ct)*nut2/Prt));
    volScalarField kByCp2("kByCp2", alpha2*(k2/Cp2/rho2 + nut2/Prt));

    fvScalarMatrix T1Eqn
    (
        fvm::ddt(alpha1, T1)
      + fvm::div(alphaPhi1, T1)
      - fvm::laplacian(kByCp1, T1)
     ==
        heatTransferCoeff*T2/Cp1/rho1
      - fvm::Sp(heatTransferCoeff/Cp1/rho1, T1)
      + alpha1*Dp1Dt/Cp1/rho1
    );

    fvScalarMatrix T2Eqn
    (
        fvm::ddt(alpha2, T2)
      + fvm::div(alphaPhi2, T2)
      - fvm::laplacian(kByCp2, T2)
     ==
        heatTransferCoeff*T1/Cp2/rho2
      - fvm::Sp(heatTransferCoeff/Cp2/rho2, T2)
      + alpha2*Dp2Dt/Cp2/rho2
    );

    T1Eqn.relax();
    T1Eqn.solve();

    T2Eqn.relax();
    T2Eqn.solve();

    // Update compressibilities
    psi1 = 1.0/(R1*T1);
    psi2 = 1.0/(R2*T2);
}
