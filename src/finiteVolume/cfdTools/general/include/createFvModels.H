Foam::fvModels& fvModels(Foam::fvModels::New(mesh));

if (!fvModels.PtrListDictionary<fvModel>::size())
{
    Info << "No fvModels present" << endl;
}
