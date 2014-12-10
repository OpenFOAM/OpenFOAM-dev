        /*
        volTensorField gradU(fvc::grad(U));
        volSymmTensorField D(symm(fvc::grad(U)));
        volTensorField Dprim(symm(fvc::grad(U - UMean)));

        volScalarField prod(-((U - UMean)*(U - UMean)) && D);
        */

        /*
        volScalarField txx
        (
            IOobject
            (
                "txx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, -1, 0, 0)
        );
        txx =sqrt(Txx - (UMeanx*UMeanx));
        txx.write();

        volScalarField tyy
        (
            IOobject
            (
                "tyy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, -1, 0, 0)
        );
        tyy = sqrt(Tyy - (UMeany*UMeany));
        tyy.write();

        volScalarField tzz
        (
            IOobject
            (
                "tzz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, -1, 0, 0)
        );
        tzz = sqrt(Tzz - (UMeanz*UMeanz));
        tzz.write();

        volScalarField txy
        (
            IOobject
            (
                "txy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, 0, 0)
        );
        txy = Txy - (UMeanx*UMeany);
        txy.write();
*/
