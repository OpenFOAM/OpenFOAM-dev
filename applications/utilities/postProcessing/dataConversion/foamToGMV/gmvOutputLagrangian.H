gmvFile << "tracers " << particles.size() << nl;

{
    pointField positions(particles.size());

    label particlei = 0;
    forAllConstIter(Cloud<passiveParticle>, particles, iter)
    {
        positions[particlei ++] = iter().position(mesh);
    }

    for (i = 0; i < pTraits<point>::nComponents; i ++)
    {
        forAll(positions, particlei)
        {
            gmvFile << component(positions[particlei], i) << ' ';
        }
        gmvFile << nl;
    }
}

forAll(lagrangianScalarNames, i)
{
    const word& name = lagrangianScalarNames[i];

    IOField<scalar> fld
    (
        IOobject
        (
            name,
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (fld.size())
    {
        gmvFile << name  << nl;

        forAll(fld, n)
        {
            gmvFile << fld[n] << token::SPACE;
        }
        gmvFile << nl;
    }
}

forAll(lagrangianVectorNames, i)
{
    const word& name = lagrangianVectorNames[i];

    IOField<vector> fld
    (
        IOobject
        (
            name,
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (fld.size())
    {
        gmvFile << name + "x" << nl;

        forAll(fld, n)
        {
            gmvFile << fld[n].x() << token::SPACE;
        }
        gmvFile << nl;

        gmvFile << name + "y" << nl;

        forAll(fld, n)
        {
            gmvFile << fld[n].y() << token::SPACE;
        }
        gmvFile << nl;

        gmvFile << name + "z" << nl;

        forAll(fld, n)
        {
            gmvFile << fld[n].z() << token::SPACE;
        }
        gmvFile << nl;
    }

}


gmvFile << "endtrace"<< nl;
