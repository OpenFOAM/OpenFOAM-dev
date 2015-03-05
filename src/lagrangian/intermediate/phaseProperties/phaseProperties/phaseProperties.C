/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "phaseProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::phaseProperties::phaseType,
        4
    >::names[] =
    {
        "gas",
        "liquid",
        "solid",
        "unknown"
    };
}

const Foam::NamedEnum<Foam::phaseProperties::phaseType, 4>
    Foam::phaseProperties::phaseTypeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::phaseProperties::reorder(const wordList& specieNames)
{
    // ***HGW Unfortunately in the current implementation it is assumed that
    // if no species are specified the phase is not present and this MUST
    // be checked at the point of use.  This needs a rewrite.
    if (!names_.size())
    {
        return;
    }

    // Store the current sames and mass-fractions
    List<word> names0(names_);
    scalarField Y0(Y_);

    // Update the specie names to those given
    names_ = specieNames;

    // Re-size mass-fractions if necessary, initialize to 0
    if (names_.size() != names0.size())
    {
        Y_.setSize(names_.size());
        Y_ = 0;
    }

    // Set the mass-fraction for each specie in the list to the corresponding
    // value in the original list
    forAll(names0, i)
    {
        bool found = false;
        forAll(names_, j)
        {
            if (names_[j] == names0[i])
            {
                Y_[j] = Y0[i];
                found = true;
                break;
            }
        }

        if (!found)
        {
            FatalErrorIn
            (
                "void phaseProperties::reorder(const wordList&)"
            )   << "Could not find specie " << names0[i]
                << " in list " <<  names_
                << " for phase " << phaseTypeNames[phase_]
                << exit(FatalError);
        }
    }
}


void Foam::phaseProperties::setCarrierIds
(
    const wordList& carrierNames
)
{
    carrierIds_ = -1;

    forAll(names_, i)
    {
        forAll(carrierNames, j)
        {
            if (carrierNames[j] == names_[i])
            {
                carrierIds_[i] = j;
                break;
            }
        }
        if (carrierIds_[i] == -1)
        {
            FatalErrorIn
            (
                "void phaseProperties::setCarrierIds"
                "(const wordList& carrierNames)"
            )   << "Could not find carrier specie " << names_[i]
                << " in species list" <<  nl
                << "Available species are: " << nl << carrierNames << nl
                << exit(FatalError);
        }
    }
}


void Foam::phaseProperties::checkTotalMassFraction() const
{
    scalar total = 0.0;
    forAll(Y_, cmptI)
    {
        total += Y_[cmptI];
    }

    if (Y_.size() != 0 && mag(total - 1.0) > SMALL)
    {
        FatalErrorIn
        (
            "void phaseProperties::checkTotalMassFraction() const"
        )   << "Component fractions must total to unity for phase "
            << phaseTypeNames[phase_] << nl
            << "Components: " << nl << names_ << nl
            << exit(FatalError);
    }
}


Foam::word Foam::phaseProperties::phaseToStateLabel(const phaseType pt) const
{
    word state = "(unknown)";
    switch (pt)
    {
        case GAS:
        {
            state = "(g)";
            break;
        }
        case LIQUID:
        {
            state = "(l)";
            break;
        }
        case SOLID:
        {
            state = "(s)";
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "phaseProperties::phaseToStateLabel(phaseType pt)"
            )   << "Invalid phase: " << phaseTypeNames[pt] << nl
                << "    phase must be gas, liquid or solid" << nl
                << exit(FatalError);
        }
    }

    return state;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseProperties::phaseProperties()
:
    phase_(UNKNOWN),
    stateLabel_("(unknown)"),
    names_(0),
    Y_(0),
    carrierIds_(0)
{}


Foam::phaseProperties::phaseProperties(const phaseProperties& pp)
:
    phase_(pp.phase_),
    stateLabel_(pp.stateLabel_),
    names_(pp.names_),
    Y_(pp.Y_),
    carrierIds_(pp.carrierIds_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseProperties::~phaseProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseProperties::reorder
(
    const wordList& gasNames,
    const wordList& liquidNames,
    const wordList& solidNames
)
{
    // Determine the addressing to map between components listed in the phase
    // with those given in the (main) thermo properties
    switch (phase_)
    {
        case GAS:
        {
            reorder(gasNames);
            forAll(carrierIds_, i)
            {
                carrierIds_[i] = i;
            }
            break;
        }
        case LIQUID:
        {
            reorder(liquidNames);
            setCarrierIds(gasNames);
            break;
        }
        case SOLID:
        {
            reorder(solidNames);
            WarningIn
            (
                "phaseProperties::reorder"
                "("
                    "const wordList& gasNames, "
                    "const wordList& liquidNames, "
                    "const wordList& solidNames"
                ")"
            )   << "Assuming no mapping between solid and carrier species"
                << endl;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "phaseProperties::reorder"
                "("
                    "const wordList& gasNames, "
                    "const wordList& liquidNames, "
                    "const wordList& solidNames"
                ")"
            )   << "Invalid phase: " << phaseTypeNames[phase_] << nl
                << "    phase must be gas, liquid or solid" << nl
                << exit(FatalError);
        }
    }
}


Foam::phaseProperties::phaseType Foam::phaseProperties::phase() const
{
    return phase_;
}


const Foam::word& Foam::phaseProperties::stateLabel() const
{
    return stateLabel_;
}


Foam::word Foam::phaseProperties::phaseTypeName() const
{
    return phaseTypeNames[phase_];
}


const Foam::List<Foam::word>& Foam::phaseProperties::names() const
{
    return names_;
}


const Foam::word& Foam::phaseProperties::name(const label cmptI) const
{
    if (cmptI >= names_.size())
    {
        FatalErrorIn
        (
            "const word& phaseProperties::name(const label) const"
        )   << "Requested component " << cmptI << "out of range" << nl
            << "Available phase components:" << nl << names_ << nl
            << exit(FatalError);
    }

    return names_[cmptI];
}


const Foam::scalarField& Foam::phaseProperties::Y() const
{
    return Y_;
}


Foam::scalar& Foam::phaseProperties::Y(const label cmptI)
{
    if (cmptI >= Y_.size())
    {
        FatalErrorIn
        (
            "const scalar& phaseProperties::Y(const label) const"
        )   << "Requested component " << cmptI << "out of range" << nl
            << "Available phase components:" << nl << names_ << nl
            << exit(FatalError);
    }

    return Y_[cmptI];
}


const Foam::labelList& Foam::phaseProperties::carrierIds() const
{
    return carrierIds_;
}


Foam::label Foam::phaseProperties::id(const word& cmptName) const
{
    forAll(names_, cmptI)
    {
        if (names_[cmptI] == cmptName)
        {
            return cmptI;
        }
    }

    return -1;
}


// ************************************************************************* //
