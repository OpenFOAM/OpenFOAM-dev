if (!mesh.foundObject<IOdictionary>("radiationProperties"))
{
    IOobject radiationProperties
    (
        "radiationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (radiationProperties.typeHeaderOk<IOdictionary>(false))
    {
        const word modelType
        (
            IOdictionary(radiationProperties).lookup("radiationModel")
        );

        if (modelType != "none")
        {
            FatalErrorInFunction
                << "Radiation model " << modelType
                << " selected but not enabled in fvOptions" << nl
                << "To enable radiation add " << nl << nl
                << "radiation" << nl
                << "{" << nl
                << "    type    radiation;" << nl
                << "    libs    (\"libradiationModels.so\");" << nl
                << "}" << nl << nl
                << "to " << radiationProperties.localPath()/"fvOptions"
                << exit(FatalError);
        }
    }
}
