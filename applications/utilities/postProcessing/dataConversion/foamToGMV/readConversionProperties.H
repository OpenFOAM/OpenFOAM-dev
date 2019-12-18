IOdictionary conversionProperties
(
    IOobject
    (
        "conversionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


scalar startTime
(
    conversionProperties.lookup<scalar>("startTime")
);

word vComp
(
    conversionProperties.lookup("vector")
);

word format
(
    conversionProperties.lookup("format")
);

if ((format != "ascii") && (format != "ieeei4r8"))
{
    FatalError
        << "format type: " << format << " unknown."
        << " Valid options are: ascii ieeei4r8"
        << abort(FatalError);

}

word cells
(
    conversionProperties.lookup("cells")
);

if
(
    (cells != "hex")
)
{
    FatalError
        << "cells type: " << cells << " unknown."
        << abort(FatalError);
}
