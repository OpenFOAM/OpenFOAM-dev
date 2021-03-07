Foam::fvConstraints& fvConstraints(Foam::fvConstraints::New(mesh));

if (!fvConstraints.PtrListDictionary<fvConstraint>::size())
{
    Info << "No fvConstraints present" << endl;
}
