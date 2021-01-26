/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "specieReactionRates.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(specieReactionRates, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        specieReactionRates,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::specieReactionRates::writeFileHeader(const label i)
{
    writeHeader(file(), "Specie reaction rates");
    volRegion::writeFileHeader(*this, file());
    writeHeaderValue(file(), "nSpecie", chemistryModel_.nSpecie());
    writeHeaderValue(file(), "nReaction", chemistryModel_.nReaction());

    writeCommented(file(), "Time");
    writeTabbed(file(), "Reaction");

    const wordList& speciesNames =
        chemistryModel_.thermo().composition().species();

    forAll (speciesNames, si)
    {
        writeTabbed(file(), speciesNames[si]);
    }

    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::specieReactionRates::specieReactionRates
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    logFiles(obr_, name),
    chemistryModel_
    (
        fvMeshFunctionObject::mesh_.lookupObject<basicChemistryModel>
        (
            "chemistryProperties"
        )
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::specieReactionRates::~specieReactionRates()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::specieReactionRates::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    resetName("specieReactionRates");

    return true;
}


bool Foam::functionObjects::specieReactionRates::execute()
{
    return true;
}


bool Foam::functionObjects::specieReactionRates::write()
{
    logFiles::write();

    const label nSpecie = chemistryModel_.nSpecie();
    const label nReaction = chemistryModel_.nReaction();

    // Region volume
    const scalar V = this->V();

    for (label ri=0; ri<nReaction; ri++)
    {
        if (Pstream::master())
        {
            writeTime(file());
            file() << token::TAB << ri;
        }

        for (label si=0; si<nSpecie; si++)
        {
            volScalarField::Internal RR
            (
                chemistryModel_.calculateRR(ri, si)
            );

            scalar sumVRRi = 0;

            if (isNull(cellIDs()))
            {
                sumVRRi = fvc::domainIntegrate(RR).value();
            }
            else
            {
                sumVRRi = gSum
                (
                    scalarField(fvMeshFunctionObject::mesh_.V()*RR, cellIDs())
                );
            }

            if (Pstream::master())
            {
                file() << token::TAB << sumVRRi/V;
            }
        }

        if (Pstream::master())
        {
            file() << nl;
        }
    }

    if (Pstream::master())
    {
        file() << nl << endl;
    }

    return true;
}


// ************************************************************************* //
