/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "processorField.H"
#include "dictionary.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(processorField, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorField::processorField
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        read(dict);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* procFieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "processorID",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless, 0.0)
            )
        );

        mesh.objectRegistry::store(procFieldPtr);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "processorField::processorField"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorField::~processorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorField::read(const dictionary& dict)
{
    // do nothing
}


void Foam::processorField::execute()
{
    if (active_)
    {
        const volScalarField& procField =
            obr_.lookupObject<volScalarField>("processorID");

        const_cast<volScalarField&>(procField) ==
            dimensionedScalar("procI", dimless, Pstream::myProcNo());
    }
}


void Foam::processorField::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::processorField::timeSet()
{
    // Do nothing
}


void Foam::processorField::write()
{
    if (active_)
    {
        const volScalarField& procField =
            obr_.lookupObject<volScalarField>("processorID");

        procField.write();
    }
}


// ************************************************************************* //
