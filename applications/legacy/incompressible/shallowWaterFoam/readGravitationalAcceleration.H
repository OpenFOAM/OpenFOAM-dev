Info<< "\nReading gravitationalProperties" << endl;

IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g
(
    "g",
    dimensions::acceleration,
    gravitationalProperties
);
const Switch rotating(gravitationalProperties.lookup("rotating"));
const dimensionedVector Omega =
    rotating
  ? dimensionedVector("Omega", dimensions::rate, gravitationalProperties)
  : dimensionedVector("Omega", dimensions::rate, vector::zero);
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;
