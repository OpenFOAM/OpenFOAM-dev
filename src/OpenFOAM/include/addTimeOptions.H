Foam::argList::addBoolOption
(
    "constant",
    "include the 'constant/' dir in the times list"
);

Foam::argList::addBoolOption
(
    "latestTime",
    "select the latest time"
);

Foam::argList::addBoolOption
(
    "noZero",
    "exclude the '0/' dir from the times list"
);

Foam::argList::addOption
(
    "time",
    "time",
    "specify a single time value to select"
);
