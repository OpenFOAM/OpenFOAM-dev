/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "vorticity.H"
#include "volFields.H"
#include "dictionary.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vorticity, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::vorticity::vorticity
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    UName_("U"),
    outputName_(typeName)
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volVectorField* vorticityPtr
    (
        new volVectorField
        (
            IOobject
            (
                outputName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("0", dimless/dimTime, Zero)
        )
    );

    mesh.objectRegistry::store(vorticityPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::vorticity::~vorticity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::vorticity::read(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("UName", "U");
    if (UName_ != "U")
    {
        outputName_ = typeName + "(" + UName_ + ")";
    }
}


void Foam::functionObjects::vorticity::execute()
{
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

    volVectorField& vorticity = const_cast<volVectorField&>
    (
        obr_.lookupObject<volVectorField>(outputName_)
    );

    vorticity = fvc::curl(U);
}


void Foam::functionObjects::vorticity::end()
{
    execute();
}


void Foam::functionObjects::vorticity::timeSet()
{}


void Foam::functionObjects::vorticity::write()
{
    const volVectorField& vorticity =
        obr_.lookupObject<volVectorField>(outputName_);

    Info<< type() << " " << name_ << " output:" << nl
        << "    writing field " << vorticity.name() << nl
        << endl;

    vorticity.write();
}


// ************************************************************************* //
