/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "reactionRates.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reactionRates, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        reactionRates,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::reactionRates::writeFileHeader(const label i)
{
    const label nReaction = chemistryModel_.nReaction();

    writeHeader(file(), "Reaction rates");

    fvCellSet::writeFileHeader(*this, file());

    writeHeaderValue(file(), "nReaction", nReaction);
    writeCommented(file(), "Time");

    for (label reactioni = 0; reactioni < nReaction; ++ reactioni)
    {
        writeTabbed(file(), chemistryModel_.reactionName(reactioni));
    }

    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reactionRates::reactionRates
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fvCellSet(fvMeshFunctionObject::mesh_, dict),
    logFiles(obr_, name),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null)),
    chemistryModel_
    (
        fvMeshFunctionObject::mesh_.lookupObject<basicChemistryModel>
        (
            IOobject::groupName("chemistryProperties", phaseName_)
        )
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::reactionRates::~reactionRates()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reactionRates::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    resetName("reactionRates");

    return true;
}


bool Foam::functionObjects::reactionRates::execute()
{
    return true;
}


bool Foam::functionObjects::reactionRates::write()
{
    logFiles::write();

    const label nReaction = chemistryModel_.nReaction();

    if (Pstream::master())
    {
        writeTime(file());
    }

    for (label reactioni = 0; reactioni < nReaction; reactioni ++)
    {
        const volScalarField::Internal RR
        (
            chemistryModel_.reactionRR(reactioni)
        );

        const scalar sumVRR =
            all()
          ? fvc::domainIntegrate(RR).value()
          : gSum
            (
                scalarField
                (
                    fvMeshFunctionObject::mesh_.V()*RR,
                    cells()
                )
            );

        if (Pstream::master())
        {
            file() << token::TAB << sumVRR/V();
        }
    }

    if (Pstream::master())
    {
        file() << endl;
    }

    return true;
}


// ************************************************************************* //
