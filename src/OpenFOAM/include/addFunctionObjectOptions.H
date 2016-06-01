#include "addDictOption.H"
Foam::argList::addOption
(
    "field",
    "name",
    "Specify the name of the field to be processed, e.g. U"
);
Foam::argList::addOption
(
    "fields",
    "list",
    "Specify a list of fields to be processed, e.g. '(U T p)' - "
    "regular expressions not currently supported"
);
Foam::argList::addOption
(
    "func",
    "name",
    "Specify the name of the functionObject to execute, e.g. Q"
);
Foam::argList::addOption
(
    "funcs",
    "list",
    "Specify the names of the functionObjects to execute, e.g. '(Q div(U))'"
);
Foam::argList::addBoolOption
(
    "list",
    "List the available configured functionObjects"
);
