if (!mesh.foundObject<IOdictionary>("radiationProperties"))
{
    typeIOobject<IOdictionary> radiationProperties
    (
        "radiationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (radiationProperties.headerOk())
    {
        const word modelType
        (
            IOdictionary(radiationProperties).lookup("radiationModel")
        );

        if (modelType != "none")
        {
            FatalErrorInFunction
                << "Radiation model " << modelType
                << " selected but not enabled in fvModels" << nl
                << "To enable radiation add " << nl << nl
                << "radiation" << nl
                << "{" << nl
                << "    type    radiation;" << nl
                << "    libs    (\"libradiationModels.so\");" << nl
                << "}" << nl << nl
                << "to " << radiationProperties.relativePath()/"fvModels"
                << exit(FatalError);
        }
    }
}
