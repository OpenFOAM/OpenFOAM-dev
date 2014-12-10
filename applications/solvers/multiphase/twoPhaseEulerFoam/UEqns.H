mrfZones.correctBoundaryVelocity(U1);
mrfZones.correctBoundaryVelocity(U2);
mrfZones.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

volScalarField dragCoeff(fluid.dragCoeff());

{
    volScalarField virtualMassCoeff(fluid.virtualMassCoeff());
    volVectorField liftForce(fluid.liftForce());
    volVectorField wallLubricationForce(fluid.wallLubricationForce());
    volVectorField turbulentDispersionForce(fluid.turbulentDispersionForce());

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(contErr1, U1)
          + mrfZones(alpha1*rho1 + virtualMassCoeff, U1)
          + phase1.turbulence().divDevRhoReff(U1)
         ==
          - liftForce
          - wallLubricationForce
          - turbulentDispersionForce
          - virtualMassCoeff
           *(
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
              - DDtU2
            )
          + fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        U1Eqn += fvm::Sp(dragCoeff, U1);
        fvOptions.constrain(U1Eqn);
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2) + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(contErr2, U2)
          + mrfZones(alpha2*rho2 + virtualMassCoeff, U2)
          + phase2.turbulence().divDevRhoReff(U2)
         ==
            liftForce
          + wallLubricationForce
          + turbulentDispersionForce
          - virtualMassCoeff
           *(
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
              - DDtU1
            )
          + fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        U2Eqn += fvm::Sp(dragCoeff, U2);
        fvOptions.constrain(U2Eqn);
    }
}
