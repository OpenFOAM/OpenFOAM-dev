Foam::fvConstraints& fvConstraints(Foam::fvConstraints::New(mesh));

if (!fvConstraints.size())
{
    Info << "No fvConstraints present" << endl;
}
