/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "chemkinReader.H"
#include <fstream>
#include "atomicWeights.H"
#include "IrreversibleReaction.H"
#include "ReversibleReaction.H"
#include "NonEquilibriumReversibleReaction.H"
#include "ArrheniusReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"
#include "FallOffReactionRate.H"
#include "ChemicallyActivatedReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "TroeFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "LandauTellerReactionRate.H"
#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"
#include "addToRunTimeSelectionTable.H"


/* * * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * */

namespace Foam
{
    addChemistryReaderType(chemkinReader, gasHThermoPhysics);
}


/* * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * */

const char* Foam::chemkinReader::reactionTypeNames[4] =
{
    "irreversible",
    "reversible",
    "nonEquilibriumReversible",
    "unknownReactionType"
};

const char* Foam::chemkinReader::reactionRateTypeNames[8] =
{
    "Arrhenius",
    "thirdBodyArrhenius",
    "unimolecularFallOff",
    "chemicallyActivatedBimolecular",
    "LandauTeller",
    "Janev",
    "powerSeries",
    "unknownReactionRateType"
};

const char* Foam::chemkinReader::fallOffFunctionNames[4] =
{
    "Lindemann",
    "Troe",
    "SRI",
    "unknownFallOffFunctionType"
};

void Foam::chemkinReader::initReactionKeywordTable()
{
    reactionKeywordTable_.insert("M", thirdBodyReactionType);
    reactionKeywordTable_.insert("LOW", unimolecularFallOffReactionType);
    reactionKeywordTable_.insert
    (
        "HIGH",
        chemicallyActivatedBimolecularReactionType
    );
    reactionKeywordTable_.insert("TROE", TroeReactionType);
    reactionKeywordTable_.insert("SRI", SRIReactionType);
    reactionKeywordTable_.insert("LT", LandauTellerReactionType);
    reactionKeywordTable_.insert("RLT", reverseLandauTellerReactionType);
    reactionKeywordTable_.insert("JAN", JanevReactionType);
    reactionKeywordTable_.insert("FIT1", powerSeriesReactionRateType);
    reactionKeywordTable_.insert("HV", radiationActivatedReactionType);
    reactionKeywordTable_.insert("TDEP", speciesTempReactionType);
    reactionKeywordTable_.insert("EXCI", energyLossReactionType);
    reactionKeywordTable_.insert("MOME", plasmaMomentumTransfer);
    reactionKeywordTable_.insert("XSMI", collisionCrossSection);
    reactionKeywordTable_.insert("REV", nonEquilibriumReversibleReactionType);
    reactionKeywordTable_.insert("DUPLICATE", duplicateReactionType);
    reactionKeywordTable_.insert("DUP", duplicateReactionType);
    reactionKeywordTable_.insert("FORD", speciesOrderForward);
    reactionKeywordTable_.insert("RORD", speciesOrderReverse);
    reactionKeywordTable_.insert("UNITS", UnitsOfReaction);
    reactionKeywordTable_.insert("END", end);
}


Foam::scalar Foam::chemkinReader::molecularWeight
(
    const List<specieElement>& specieComposition
) const
{
    scalar molWt = 0.0;

    forAll(specieComposition, i)
    {
        label nAtoms = specieComposition[i].nAtoms;
        const word& elementName = specieComposition[i].elementName;

        if (isotopeAtomicWts_.found(elementName))
        {
            molWt += nAtoms*isotopeAtomicWts_[elementName];
        }
        else if (atomicWeights.found(elementName))
        {
            molWt += nAtoms*atomicWeights[elementName];
        }
        else
        {
            FatalErrorIn("chemkinReader::lex()")
                << "Unknown element " << elementName
                << " on line " << lineNo_-1 << nl
                << "    specieComposition: " << specieComposition
                << exit(FatalError);
        }
    }

    return molWt;
}


void Foam::chemkinReader::checkCoeffs
(
    const scalarList& reactionCoeffs,
    const char* reactionRateName,
    const label nCoeffs
) const
{
    if (reactionCoeffs.size() != nCoeffs)
    {
        FatalErrorIn("chemkinReader::checkCoeffs")
            << "Wrong number of coefficients for the " << reactionRateName
            << " rate expression on line "
            << lineNo_-1 << ", should be "
            << nCoeffs << " but " << reactionCoeffs.size() << " supplied." << nl
            << "Coefficients are "
            << reactionCoeffs << nl
            << exit(FatalError);
    }
}

template<class ReactionRateType>
void Foam::chemkinReader::addReactionType
(
    const reactionType rType,
    DynamicList<gasHReaction::specieCoeffs>& lhs,
    DynamicList<gasHReaction::specieCoeffs>& rhs,
    const ReactionRateType& rr
)
{
    switch (rType)
    {
        case irreversible:
        {
            reactions_.append
            (
                new IrreversibleReaction
                <Reaction, gasHThermoPhysics, ReactionRateType>
                (
                    Reaction<gasHThermoPhysics>
                    (
                        speciesTable_,
                        lhs.shrink(),
                        rhs.shrink(),
                        speciesThermo_
                    ),
                    rr
                )
            );
        }
        break;

        case reversible:
        {
            reactions_.append
            (
                new ReversibleReaction
                <Reaction, gasHThermoPhysics, ReactionRateType>
                (
                    Reaction<gasHThermoPhysics>
                    (
                        speciesTable_,
                        lhs.shrink(),
                        rhs.shrink(),
                        speciesThermo_
                    ),
                    rr
                )
            );
        }
        break;

        default:

            if (rType < 3)
            {
                FatalErrorIn("chemkinReader::addReactionType")
                    << "Reaction type " << reactionTypeNames[rType]
                    << " on line " << lineNo_-1
                    << " not handled by this function"
                    << exit(FatalError);
            }
            else
            {
                FatalErrorIn("chemkinReader::addReactionType")
                    << "Unknown reaction type " << rType
                    << " on line " << lineNo_-1
                    << exit(FatalError);
            }
    }
}

template<template<class, class> class PressureDependencyType>
void Foam::chemkinReader::addPressureDependentReaction
(
    const reactionType rType,
    const fallOffFunctionType fofType,
    DynamicList<gasHReaction::specieCoeffs>& lhs,
    DynamicList<gasHReaction::specieCoeffs>& rhs,
    const scalarList& efficiencies,
    const scalarList& k0Coeffs,
    const scalarList& kInfCoeffs,
    const HashTable<scalarList>& reactionCoeffsTable,
    const scalar Afactor0,
    const scalar AfactorInf,
    const scalar RR
)
{
    checkCoeffs(k0Coeffs, "k0", 3);
    checkCoeffs(kInfCoeffs, "kInf", 3);

    switch (fofType)
    {
        case Lindemann:
        {
            addReactionType
            (
                rType,
                lhs, rhs,
                PressureDependencyType
                    <ArrheniusReactionRate, LindemannFallOffFunction>
                (
                    ArrheniusReactionRate
                    (
                        Afactor0*k0Coeffs[0],
                        k0Coeffs[1],
                        k0Coeffs[2]/RR
                    ),
                    ArrheniusReactionRate
                    (
                        AfactorInf*kInfCoeffs[0],
                        kInfCoeffs[1],
                        kInfCoeffs[2]/RR
                    ),
                    LindemannFallOffFunction(),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
            break;
        }
        case Troe:
        {
            scalarList TroeCoeffs
            (
                reactionCoeffsTable[fallOffFunctionNames[fofType]]
            );

            if (TroeCoeffs.size() != 4 && TroeCoeffs.size() != 3)
            {
                FatalErrorIn("chemkinReader::addPressureDependentReaction")
                    << "Wrong number of coefficients for Troe rate expression"
                       " on line " << lineNo_-1 << ", should be 3 or 4 but "
                    << TroeCoeffs.size() << " supplied." << nl
                    << "Coefficients are "
                    << TroeCoeffs << nl
                    << exit(FatalError);
            }

            if (TroeCoeffs.size() == 3)
            {
                TroeCoeffs.setSize(4);
                TroeCoeffs[3] = GREAT;
            }

            addReactionType
            (
                rType,
                lhs, rhs,
                PressureDependencyType
                    <ArrheniusReactionRate, TroeFallOffFunction>
                (
                    ArrheniusReactionRate
                    (
                        Afactor0*k0Coeffs[0],
                        k0Coeffs[1],
                        k0Coeffs[2]/RR
                    ),
                    ArrheniusReactionRate
                    (
                        AfactorInf*kInfCoeffs[0],
                        kInfCoeffs[1],
                        kInfCoeffs[2]/RR
                    ),
                    TroeFallOffFunction
                    (
                        TroeCoeffs[0],
                        TroeCoeffs[1],
                        TroeCoeffs[2],
                        TroeCoeffs[3]
                    ),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
            break;
        }
        case SRI:
        {
            scalarList SRICoeffs
            (
                reactionCoeffsTable[fallOffFunctionNames[fofType]]
            );

            if (SRICoeffs.size() != 5 && SRICoeffs.size() != 3)
            {
                FatalErrorIn("chemkinReader::addPressureDependentReaction")
                    << "Wrong number of coefficients for SRI rate expression"
                       " on line " << lineNo_-1 << ", should be 3 or 5 but "
                    << SRICoeffs.size() << " supplied." << nl
                    << "Coefficients are "
                    << SRICoeffs << nl
                    << exit(FatalError);
            }

            if (SRICoeffs.size() == 3)
            {
                SRICoeffs.setSize(5);
                SRICoeffs[3] = 1.0;
                SRICoeffs[4] = 0.0;
            }

            addReactionType
            (
                rType,
                lhs, rhs,
                PressureDependencyType
                    <ArrheniusReactionRate, SRIFallOffFunction>
                (
                    ArrheniusReactionRate
                    (
                        Afactor0*k0Coeffs[0],
                        k0Coeffs[1],
                        k0Coeffs[2]/RR
                    ),
                    ArrheniusReactionRate
                    (
                        AfactorInf*kInfCoeffs[0],
                        kInfCoeffs[1],
                        kInfCoeffs[2]/RR
                    ),
                    SRIFallOffFunction
                    (
                        SRICoeffs[0],
                        SRICoeffs[1],
                        SRICoeffs[2],
                        SRICoeffs[3],
                        SRICoeffs[4]
                    ),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
            break;
        }
        default:
        {
            FatalErrorIn("chemkinReader::addPressureDependentReaction")
                << "Fall-off function type "
                << fallOffFunctionNames[fofType]
                << " on line " << lineNo_-1
                << " not implemented"
                << exit(FatalError);
        }
    }
}


void Foam::chemkinReader::addReaction
(
    DynamicList<gasHReaction::specieCoeffs>& lhs,
    DynamicList<gasHReaction::specieCoeffs>& rhs,
    const scalarList& efficiencies,
    const reactionType rType,
    const reactionRateType rrType,
    const fallOffFunctionType fofType,
    const scalarList& ArrheniusCoeffs,
    HashTable<scalarList>& reactionCoeffsTable,
    const scalar RR
)
{
    checkCoeffs(ArrheniusCoeffs, "Arrhenius", 3);

    scalarList nAtoms(elementNames_.size(), 0.0);

    forAll(lhs, i)
    {
        const List<specieElement>& specieComposition =
            specieComposition_[speciesTable_[lhs[i].index]];

        forAll(specieComposition, j)
        {
            label elementi = elementIndices_[specieComposition[j].elementName];
            nAtoms[elementi] += lhs[i].stoichCoeff*specieComposition[j].nAtoms;
        }
    }

    forAll(rhs, i)
    {
        const List<specieElement>& specieComposition =
            specieComposition_[speciesTable_[rhs[i].index]];

        forAll(specieComposition, j)
        {
            label elementi = elementIndices_[specieComposition[j].elementName];
            nAtoms[elementi] -= rhs[i].stoichCoeff*specieComposition[j].nAtoms;
        }
    }


    // Calculate the unit conversion factor for the A coefficient
    // for the change from mol/cm^3 to kmol/m^3 concentraction units
    const scalar concFactor = 0.001;
    scalar sumExp = 0.0;
    forAll(lhs, i)
    {
        sumExp += lhs[i].exponent;
    }
    scalar Afactor = pow(concFactor, sumExp - 1.0);

    scalar AfactorRev = Afactor;

    if (rType == nonEquilibriumReversible)
    {
        sumExp = 0.0;
        forAll(rhs, i)
        {
            sumExp += rhs[i].exponent;
        }
        AfactorRev = pow(concFactor, sumExp - 1.0);
    }

    switch (rrType)
    {
        case Arrhenius:
        {
            if (rType == nonEquilibriumReversible)
            {
                const scalarList& reverseArrheniusCoeffs =
                    reactionCoeffsTable[reactionTypeNames[rType]];

                checkCoeffs(reverseArrheniusCoeffs, "reverse Arrhenius", 3);

                reactions_.append
                (
                    new NonEquilibriumReversibleReaction
                        <Reaction, gasHThermoPhysics, ArrheniusReactionRate>
                    (
                        Reaction<gasHThermoPhysics>
                        (
                            speciesTable_,
                            lhs.shrink(),
                            rhs.shrink(),
                            speciesThermo_
                        ),
                        ArrheniusReactionRate
                        (
                            Afactor*ArrheniusCoeffs[0],
                            ArrheniusCoeffs[1],
                            ArrheniusCoeffs[2]/RR
                        ),
                        ArrheniusReactionRate
                        (
                            AfactorRev*reverseArrheniusCoeffs[0],
                            reverseArrheniusCoeffs[1],
                            reverseArrheniusCoeffs[2]/RR
                        )
                    )
                );
            }
            else
            {
                addReactionType
                (
                    rType,
                    lhs, rhs,
                    ArrheniusReactionRate
                    (
                        Afactor*ArrheniusCoeffs[0],
                        ArrheniusCoeffs[1],
                        ArrheniusCoeffs[2]/RR
                    )
                );
            }
            break;
        }
        case thirdBodyArrhenius:
        {
            if (rType == nonEquilibriumReversible)
            {
                const scalarList& reverseArrheniusCoeffs =
                    reactionCoeffsTable[reactionTypeNames[rType]];

                checkCoeffs(reverseArrheniusCoeffs, "reverse Arrhenius", 3);

                reactions_.append
                (
                    new NonEquilibriumReversibleReaction
                    <
                        Reaction,
                        gasHThermoPhysics,
                        thirdBodyArrheniusReactionRate
                    >
                    (
                        Reaction<gasHThermoPhysics>
                        (
                            speciesTable_,
                            lhs.shrink(),
                            rhs.shrink(),
                            speciesThermo_
                        ),
                        thirdBodyArrheniusReactionRate
                        (
                            Afactor*concFactor*ArrheniusCoeffs[0],
                            ArrheniusCoeffs[1],
                            ArrheniusCoeffs[2]/RR,
                            thirdBodyEfficiencies(speciesTable_, efficiencies)
                        ),
                        thirdBodyArrheniusReactionRate
                        (
                            AfactorRev*concFactor*reverseArrheniusCoeffs[0],
                            reverseArrheniusCoeffs[1],
                            reverseArrheniusCoeffs[2]/RR,
                            thirdBodyEfficiencies(speciesTable_, efficiencies)
                        )
                    )
                );
            }
            else
            {
                addReactionType
                (
                    rType,
                    lhs, rhs,
                    thirdBodyArrheniusReactionRate
                    (
                        Afactor*concFactor*ArrheniusCoeffs[0],
                        ArrheniusCoeffs[1],
                        ArrheniusCoeffs[2]/RR,
                        thirdBodyEfficiencies(speciesTable_, efficiencies)
                    )
                );
            }
            break;
        }
        case unimolecularFallOff:
        {
            addPressureDependentReaction<FallOffReactionRate>
            (
                rType,
                fofType,
                lhs,
                rhs,
                efficiencies,
                reactionCoeffsTable[reactionRateTypeNames[rrType]],
                ArrheniusCoeffs,
                reactionCoeffsTable,
                concFactor*Afactor,
                Afactor,
                RR
            );
            break;
        }
        case chemicallyActivatedBimolecular:
        {
            addPressureDependentReaction<ChemicallyActivatedReactionRate>
            (
                rType,
                fofType,
                lhs,
                rhs,
                efficiencies,
                ArrheniusCoeffs,
                reactionCoeffsTable[reactionRateTypeNames[rrType]],
                reactionCoeffsTable,
                Afactor,
                Afactor/concFactor,
                RR
            );
            break;
        }
        case LandauTeller:
        {
            const scalarList& LandauTellerCoeffs =
                reactionCoeffsTable[reactionRateTypeNames[rrType]];
            checkCoeffs(LandauTellerCoeffs, "Landau-Teller", 2);

            if (rType == nonEquilibriumReversible)
            {
                const scalarList& reverseArrheniusCoeffs =
                    reactionCoeffsTable[reactionTypeNames[rType]];
                checkCoeffs(reverseArrheniusCoeffs, "reverse Arrhenius", 3);

                const scalarList& reverseLandauTellerCoeffs =
                    reactionCoeffsTable
                    [
                        word(reactionTypeNames[rType])
                      + reactionRateTypeNames[rrType]
                    ];
                checkCoeffs(LandauTellerCoeffs, "reverse Landau-Teller", 2);

                reactions_.append
                (
                    new NonEquilibriumReversibleReaction
                        <Reaction, gasHThermoPhysics, LandauTellerReactionRate>
                    (
                        Reaction<gasHThermoPhysics>
                        (
                            speciesTable_,
                            lhs.shrink(),
                            rhs.shrink(),
                            speciesThermo_
                        ),
                        LandauTellerReactionRate
                        (
                            Afactor*ArrheniusCoeffs[0],
                            ArrheniusCoeffs[1],
                            ArrheniusCoeffs[2]/RR,
                            LandauTellerCoeffs[0],
                            LandauTellerCoeffs[1]
                        ),
                        LandauTellerReactionRate
                        (
                            AfactorRev*reverseArrheniusCoeffs[0],
                            reverseArrheniusCoeffs[1],
                            reverseArrheniusCoeffs[2]/RR,
                            reverseLandauTellerCoeffs[0],
                            reverseLandauTellerCoeffs[1]
                        )
                    )
                );
            }
            else
            {
                addReactionType
                (
                    rType,
                    lhs, rhs,
                    LandauTellerReactionRate
                    (
                        Afactor*ArrheniusCoeffs[0],
                        ArrheniusCoeffs[1],
                        ArrheniusCoeffs[2]/RR,
                        LandauTellerCoeffs[0],
                        LandauTellerCoeffs[1]
                    )
                );
            }
            break;
        }
        case Janev:
        {
            const scalarList& JanevCoeffs =
                reactionCoeffsTable[reactionRateTypeNames[rrType]];

            checkCoeffs(JanevCoeffs, "Janev", 9);

            addReactionType
            (
                rType,
                lhs, rhs,
                JanevReactionRate
                (
                    Afactor*ArrheniusCoeffs[0],
                    ArrheniusCoeffs[1],
                    ArrheniusCoeffs[2]/RR,
                    FixedList<scalar, 9>(JanevCoeffs)
                )
            );
            break;
        }
        case powerSeries:
        {
            const scalarList& powerSeriesCoeffs =
                reactionCoeffsTable[reactionRateTypeNames[rrType]];

            checkCoeffs(powerSeriesCoeffs, "power-series", 4);

            addReactionType
            (
                rType,
                lhs, rhs,
                powerSeriesReactionRate
                (
                    Afactor*ArrheniusCoeffs[0],
                    ArrheniusCoeffs[1],
                    ArrheniusCoeffs[2]/RR,
                    FixedList<scalar, 4>(powerSeriesCoeffs)
                )
            );
            break;
        }
        case unknownReactionRateType:
        {
            FatalErrorIn("chemkinReader::addReaction")
                << "Internal error on line " << lineNo_-1
                << ": reaction rate type has not been set"
                << exit(FatalError);
            break;
        }
        default:
        {
            FatalErrorIn("chemkinReader::addReaction")
                << "Reaction rate type " << reactionRateTypeNames[rrType]
                << " on line " << lineNo_-1
                << " not implemented"
                << exit(FatalError);
        }
    }


    forAll(nAtoms, i)
    {
        if (mag(nAtoms[i]) > SMALL)
        {
            FatalErrorIn("chemkinReader::addReaction")
                << "Elemental imbalance in " << elementNames_[i]
                << " in reaction" << nl
                << reactions_.last() << nl
                << " on line " << lineNo_-1
                << exit(FatalError);
        }
    }

    lhs.clear();
    rhs.clear();
    reactionCoeffsTable.clear();
}


void Foam::chemkinReader::read
(
    const fileName& CHEMKINFileName,
    const fileName& thermoFileName
)
{
    if (thermoFileName != fileName::null)
    {
        std::ifstream thermoStream(thermoFileName.c_str());

        if (!thermoStream)
        {
            FatalErrorIn
            (
                "chemkin::chemkin(const fileName& CHEMKINFileName, "
                "const fileName& thermoFileName)"
            )   << "file " << thermoFileName << " not found"
                << exit(FatalError);
        }

        yy_buffer_state* bufferPtr(yy_create_buffer(&thermoStream, yyBufSize));
        yy_switch_to_buffer(bufferPtr);

        while (lex() != 0)
        {}

        yy_delete_buffer(bufferPtr);

        lineNo_ = 1;
    }

    std::ifstream CHEMKINStream(CHEMKINFileName.c_str());

    if (!CHEMKINStream)
    {
        FatalErrorIn
        (
            "chemkin::chemkin(const fileName& CHEMKINFileName, "
            "const fileName& thermoFileName)"
        )   << "file " << CHEMKINFileName << " not found"
            << exit(FatalError);
    }

    yy_buffer_state* bufferPtr(yy_create_buffer(&CHEMKINStream, yyBufSize));
    yy_switch_to_buffer(bufferPtr);

    initReactionKeywordTable();

    while (lex() != 0)
    {}

    yy_delete_buffer(bufferPtr);
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::chemkinReader::chemkinReader
(
    const fileName& CHEMKINFileName,
    speciesTable& species,
    const fileName& thermoFileName,
    const bool newFormat
)
:
    lineNo_(1),
    specieNames_(10),
    speciesTable_(species),
    reactions_(speciesTable_, speciesThermo_),
    newFormat_(newFormat)
{
    read(CHEMKINFileName, thermoFileName);
}


Foam::chemkinReader::chemkinReader
(
    const dictionary& thermoDict,
    speciesTable& species
)
:
    lineNo_(1),
    specieNames_(10),
    speciesTable_(species),
    reactions_(speciesTable_, speciesThermo_),
    newFormat_(thermoDict.lookupOrDefault("newFormat", false))
{
    if (newFormat_)
    {
        Info<< "Reading CHEMKIN thermo data in new file format" << endl;
    }

    fileName chemkinFile(fileName(thermoDict.lookup("CHEMKINFile")).expand());

    fileName thermoFile = fileName::null;

    if (thermoDict.found("CHEMKINThermoFile"))
    {
        thermoFile = fileName(thermoDict.lookup("CHEMKINThermoFile")).expand();
    }

    // allow relative file names
    fileName relPath = thermoDict.name().path();
    if (relPath.size())
    {
        if (!chemkinFile.isAbsolute())
        {
            chemkinFile = relPath/chemkinFile;
        }

        if (!thermoFile.isAbsolute())
        {
            thermoFile = relPath/thermoFile;
        }
    }

    read(chemkinFile, thermoFile);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
