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

const dimensionedVector g(gravitationalProperties.lookup("g"));
const Switch rotating(gravitationalProperties.lookup("rotating"));
const dimensionedVector Omega =
    rotating
  ? gravitationalProperties.lookup("Omega")
  : dimensionedVector("Omega", dimless/dimTime, vector::zero);
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;
