/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
#include "basicChemistryModel.H"
#include "fvcVolumeIntegrate.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
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
    const basicChemistryModel& chemistryModel =
        mesh().lookupObject<basicChemistryModel>
        (
            IOobject::groupName("chemistryProperties", phaseName_)
        );

    const label nReaction = chemistryModel.nReaction();

    writeHeader(file(), "Reaction rates");

    zone_.writeFileHeader(*this, file());

    writeHeaderValue(file(), "nReaction", nReaction);
    writeCommented(file(), "Time");

    for (label reactioni = 0; reactioni < nReaction; ++ reactioni)
    {
        writeTabbed(file(), chemistryModel.reactionName(reactioni));
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
    logFiles(obr_, name),
    zone_(fvMeshFunctionObject::mesh_, dict),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null)),
    writeFields_(false)
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

    writeFields_ = dict.lookupOrDefault<bool>("writeFields", false);

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

    const basicChemistryModel& chemistryModel =
        mesh().lookupObject<basicChemistryModel>
        (
            IOobject::groupName("chemistryProperties", phaseName_)
        );

    const label nReaction = chemistryModel.nReaction();

    if (Pstream::master())
    {
        writeTime(file());
    }

    for (label reactioni = 0; reactioni < nReaction; reactioni ++)
    {
        const volScalarField::Internal RR
        (
            chemistryModel.reactionRR(reactioni)
        );

        // Compute the average rate and write it into the log file
        const scalar sumVRR =
            zone_.all()
          ? fvc::domainIntegrate(RR).value()
          : gSum
            (
                scalarField
                (
                    fvMeshFunctionObject::mesh_.V()*RR,
                    zone_.zone()
                )
            );

        if (Pstream::master())
        {
            file() << token::TAB << sumVRR/zone_.V();
        }

        // Write the rate field, if necessary
        if (writeFields_)
        {
            RR.write();
        }
    }

    if (Pstream::master())
    {
        file() << endl;
    }

    return true;
}


void Foam::functionObjects::reactionRates::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &this->mesh())
    {
        zone_.movePoints();
    }
}


void Foam::functionObjects::reactionRates::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.topoChange(map);
    }
}


void Foam::functionObjects::reactionRates::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.mapMesh(map);
    }
}


void Foam::functionObjects::reactionRates::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.distribute(map);
    }
}


// ************************************************************************* //
