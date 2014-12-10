gmvFile << "tracers " << particles.size() << nl;
forAllConstIter(Cloud<passiveParticle>, particles, iter)
{
    gmvFile << iter().position().x() << ' ';
}
gmvFile << nl;

forAllConstIter(Cloud<passiveParticle>, particles, iter)
{
    gmvFile << iter().position().y() << ' ';
}
gmvFile << nl;

forAllConstIter(Cloud<passiveParticle>, particles, iter)
{
    gmvFile << iter().position().z() << ' ';
}
gmvFile << nl;

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
