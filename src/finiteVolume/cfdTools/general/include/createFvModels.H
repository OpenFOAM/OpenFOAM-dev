Foam::fvModels& fvModels(Foam::fvModels::New(mesh));

if (!fvModels.size())
{
    Info << "No fvModels present" << endl;
}
