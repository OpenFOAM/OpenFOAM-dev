/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
#include "volFields.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ChemistryModelType>
void Foam::functionObjects::specieReactionRates<ChemistryModelType>::
writeFileHeader
(
    const label i
)
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

template<class ChemistryModelType>
Foam::functionObjects::specieReactionRates<ChemistryModelType>::
specieReactionRates
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
        fvMeshFunctionObject::mesh_.lookupObject<ChemistryModelType>
        (
            "chemistryProperties"
        )
    )
{
    resetName("specieReactionRates");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModelType>
Foam::functionObjects::specieReactionRates<ChemistryModelType>::
~specieReactionRates()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::read
(
    const dictionary& dict
)
{
    regionFunctionObject::read(dict);

    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::execute()
{
    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::write()
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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#include "addToRunTimeSelectionTable.H"
#include "BasicChemistryModel.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

namespace Foam
{
    typedef
        functionObjects::specieReactionRates
        <
            BasicChemistryModel
            <
                psiReactionThermo
            >
        >
        psiSpecieReactionRates;

    defineTemplateTypeNameAndDebugWithName
    (
        psiSpecieReactionRates,
        "psiSpecieReactionRates",
        0
    );

    addToRunTimeSelectionTable
    (
        functionObject,
        psiSpecieReactionRates,
        dictionary
    );


    typedef
        functionObjects::specieReactionRates
        <
            BasicChemistryModel
            <
                rhoReactionThermo
            >
        >
        rhoSpecieReactionRates;

    defineTemplateTypeNameAndDebugWithName
    (
        rhoSpecieReactionRates,
        "rhoSpecieReactionRates",
        0
    );

    addToRunTimeSelectionTable
    (
        functionObject,
        rhoSpecieReactionRates,
        dictionary
    );
}


// ************************************************************************* //
