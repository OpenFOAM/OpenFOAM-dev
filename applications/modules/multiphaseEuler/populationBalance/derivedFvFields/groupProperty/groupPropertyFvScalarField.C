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

#include "groupPropertyFvScalarField.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::groupPropertyFvScalarField::groupPropertyFvScalarField
(
    const DimensionedField<scalar, volMesh>& iF
)
:
    internalField_(iF),
    popBalPtr_(nullptr),
    i_(-1)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::populationBalanceModel&
Foam::groupPropertyFvScalarField::popBal() const
{
    // Use the cached pointer if we have it
    if (popBalPtr_ != nullptr) return *popBalPtr_;

    // Look up the phase system
    const phaseSystem& fluid =
        internalField_.db().lookupObject<phaseSystem>
        (
            phaseSystem::propertiesName
        );

    // Get the diameter model for this phase and check its type
    const diameterModel& diameter =
        fluid.phases()[internalField_.group()].diameter();
    if (!isA<diameterModels::populationBalance>(diameter))
    {
        FatalErrorInFunction
            << "Phase " << internalField_.group() << " does not have a "
            << " diameter model associated with a population balance model"
            << exit(FatalError);
    }

    // Cast the diameter model and extract the population balance model
    popBalPtr_ =
        &refCast<const diameterModels::populationBalance>(diameter).popBal();

    return *popBalPtr_;
}


Foam::label Foam::groupPropertyFvScalarField::i() const
{
    // Use the cached index if we have it
    if (i_ != -1) return i_;

    // Remove the phase suffix
    const word member = internalField_.member();

    // Determine the number of trailing digits
    string::size_type memberChari = member.size();
    while (memberChari && isdigit(member[memberChari - 1]))
    {
        memberChari --;
    }

    // Check that we have something
    if (memberChari == member.size())
    {
        FatalErrorInFunction
            << "Could not determine the group index for the property field "
            << internalField_.name() << exit(FatalError);
    }

    // Convert the trailing digits to an integer
    i_ = atoi(member(memberChari, member.size()).c_str());

    // Check that the index is in bounds
    if (i_ < 0 || i_ >= popBal().nGroups())
    {
        FatalErrorInFunction
            << "Group index for property field " << internalField_.name()
            << " is out of range" << exit(FatalError);
    }

    return i_;
}


// ************************************************************************* //
