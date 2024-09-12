/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

    fvCellSet::writeFileHeader(*this, file());

    writeHeaderValue(file(), "nSpecie", chemistryModel_.nSpecie());
    writeHeaderValue(file(), "nReaction", chemistryModel_.nReaction());

    writeCommented(file(), "Time");
    writeTabbed(file(), "Reaction");

    const wordList& speciesNames =
        chemistryModel_.thermo().species();

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
    fvCellSet(fvMeshFunctionObject::mesh_, dict),
    logFiles(obr_, name),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null)),
    chemistryModel_
    (
        fvMeshFunctionObject::mesh_.lookupObject<basicChemistryModel>
        (
            IOobject::groupName("chemistryProperties", phaseName_)
        )
    ),
    writeFields_(false)
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

    writeFields_ = dict.lookupOrDefault<bool>("writeFields", false);

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

    for (label reactioni=0; reactioni<nReaction; reactioni++)
    {
        if (Pstream::master())
        {
            writeTime(file());
            file() << token::TAB << reactioni;
        }

        const PtrList<volScalarField::Internal> RR
        (
            chemistryModel_.specieReactionRR(reactioni)
        );

        // Compute the average rates and write them into the log file
        for (label speciei=0; speciei<nSpecie; speciei++)
        {
            scalar sumVRRi = 0;

            if (all())
            {
                sumVRRi = fvc::domainIntegrate(RR[speciei]).value();
            }
            else
            {
                sumVRRi =
                    gSum
                    (
                        scalarField
                        (
                            fvMeshFunctionObject::mesh_.V()*RR[speciei],
                            cells()
                        )
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

        // Write the rate fields, if necessary
        if (writeFields_)
        {
            for (label speciei=0; speciei<nSpecie; speciei++)
            {
                RR[speciei].write();
            }
        }
    }

    if (Pstream::master())
    {
        file() << endl;
    }

    return true;
}


// ************************************************************************* //
