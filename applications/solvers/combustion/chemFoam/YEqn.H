{
    forAll(Y, specieI)
    {
        volScalarField& Yi = Y[specieI];

        solve
        (
            fvm::ddt(rho, Yi) - chemistry.RR(specieI),
            mesh.solver("Yi")
        );
    }
}