#include "addDictOption.H"
Foam::argList::addOption
(
    "field",
    "name",
    "specify the name of the field to be processed, e.g. U"
);
Foam::argList::addOption
(
    "fields",
    "list",
    "specify a list of fields to be processed, e.g. '(U T p)' - "
    "regular expressions not currently supported"
);
Foam::argList::addOption
(
    "func",
    "name",
    "specify the name of the functionObject to execute, e.g. Q"
);
Foam::argList::addOption
(
    "funcs",
    "list",
    "specify the names of the functionObjects to execute, e.g. '(Q div(U))'"
);
