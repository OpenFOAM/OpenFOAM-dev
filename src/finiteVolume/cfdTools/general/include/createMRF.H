const Foam::MRFZones& MRF(Foam::MRFZones::New(mesh));

if (!MRF.size())
{
    Info << "No MRF regions present" << endl;
}
