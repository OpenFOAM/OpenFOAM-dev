#include "listOptions.H"

Foam::argList args(argc, argv);
if (!args.checkRootCase())
{
    Foam::FatalError.exit();
}

#include "listOutput.H"
