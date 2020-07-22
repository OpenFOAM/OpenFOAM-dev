/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
            FatalErrorInFunction
                << "Could not find specie " << names0[i]
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
            FatalErrorInFunction
                << "Could not find carrier specie " << names_[i]
                << " in species list" <<  nl
                << "Available species are: " << nl << carrierNames << nl
                << exit(FatalError);
        }
    }
}


void Foam::phaseProperties::checkTotalMassFraction() const
{
    scalar total = 0.0;
    forAll(Y_, speciei)
    {
        total += Y_[speciei];
    }

    if (Y_.size() != 0 && mag(total - 1.0) > small)
    {
        FatalErrorInFunction
            << "Specie fractions must total to unity for phase "
            << phaseTypeNames[phase_] << nl
            << "Species: " << nl << names_ << nl
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
            FatalErrorInFunction
                << "Invalid phase: " << phaseTypeNames[pt] << nl
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
    // Determine the addressing to map between species listed in the phase
    // with those given in the (main) thermo properties
    switch (phase_)
    {
        case GAS:
        {
            // The list of gaseous species in the mixture may be a sub-set of
            // the gaseous species in the carrier phase
            setCarrierIds(gasNames);
            break;
        }
        case LIQUID:
        {
            // Set the list of liquid species to correspond to the complete list
            // defined in the thermodynamics package.
            reorder(liquidNames);
            // Set the ids of the corresponding species in the carrier phase
            setCarrierIds(gasNames);
            break;
        }
        case SOLID:
        {
            // Set the list of solid species to correspond to the complete list
            // defined in the thermodynamics package.
            reorder(solidNames);
            // Assume there is no correspondence between the solid species and
            // the species in the carrier phase (no sublimation).
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Invalid phase: " << phaseTypeNames[phase_] << nl
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


const Foam::word& Foam::phaseProperties::name(const label speciei) const
{
    if (speciei >= names_.size())
    {
        FatalErrorInFunction
            << "Requested specie " << speciei << "out of range" << nl
            << "Available phase species:" << nl << names_ << nl
            << exit(FatalError);
    }

    return names_[speciei];
}


const Foam::scalarField& Foam::phaseProperties::Y() const
{
    return Y_;
}


Foam::scalar& Foam::phaseProperties::Y(const label speciei)
{
    if (speciei >= Y_.size())
    {
        FatalErrorInFunction
            << "Requested specie " << speciei << "out of range" << nl
            << "Available phase species:" << nl << names_ << nl
            << exit(FatalError);
    }

    return Y_[speciei];
}


const Foam::labelList& Foam::phaseProperties::carrierIds() const
{
    return carrierIds_;
}


Foam::label Foam::phaseProperties::id(const word& specieName) const
{
    forAll(names_, speciei)
    {
        if (names_[speciei] == specieName)
        {
            return speciei;
        }
    }

    return -1;
}


// ************************************************************************* //
