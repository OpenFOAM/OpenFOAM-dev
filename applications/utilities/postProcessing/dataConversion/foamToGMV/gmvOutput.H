for (label i=0; i < nTypes; i++)
{
    wordList fieldNames = objects.names(fieldTypes[i]);

    if (fieldTypes[i] == "volScalarField")
    {
        gmvFile << "variable" << nl;
    }

    forAll(fieldNames, j)
    {
        const word& fieldName = fieldNames[j];

        IOobject fieldObject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldTypes[i] == "volScalarField")
        {
            volScalarField fld(fieldObject, mesh);
            gmvFile << fieldName << " 0" << nl;
            for (label indx=0;indx<mesh.nCells();indx++)
            {
                gmvFile << fld[indx] << " ";
            }
            gmvFile << nl;
        }

        if (fieldTypes[i] == "volVectorField")
        {
            if (fieldName == vComp)
            {
                volVectorField fld(fieldObject, mesh);
                gmvFile << "velocity 0" << nl;
                for (label indx=0;indx<mesh.nCells();indx++)
                {
                    gmvFile << fld[indx].x() << " ";
                }
                for (label indx=0;indx<mesh.nCells();indx++)
                {
                    gmvFile << fld[indx].y() << " ";
                }
                for (label indx=0;indx<mesh.nCells();indx++)
                {
                    gmvFile << fld[indx].z() << " ";
                }
                gmvFile << nl;
            }
        }

        if (fieldTypes[i] == "surfaceScalarField")
        {
            // ...
        }

    }

    if (fieldTypes[i] == cloud::prefix)
    {
        IOobject lagrangianHeader
        (
            "positions",
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::NO_READ
        );

        if (lagrangianHeader.headerOk())
        {
            Cloud<passiveParticle> particles(mesh);

            IOobjectList objects(mesh, runTime.timeName(), cloud::prefix);

            wordList lagrangianScalarNames = objects.names("scalarField");
            wordList lagrangianVectorNames = objects.names("vectorField");

            if (particles.size())
            {
#               include "gmvOutputLagrangian.H"
            }
        }
    }

    if (fieldTypes[i] == "volScalarField")
    {
        gmvFile << "endvars" << nl;
    }
}
