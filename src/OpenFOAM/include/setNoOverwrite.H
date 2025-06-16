const bool overwrite = !args.optionFound("noOverwrite");

if (args.optionFound("overwrite"))
{
    Foam::Info
        << "Deprecated option 'overwrite' specified, "
           "this is now the default behaviour"
        << Foam::nl << Foam::endl;
}
