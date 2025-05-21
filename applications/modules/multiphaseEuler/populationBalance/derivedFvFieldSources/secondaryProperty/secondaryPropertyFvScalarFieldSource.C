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

#include "secondaryPropertyFvScalarFieldSource.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::secondaryPropertyFvScalarFieldSource::
secondaryPropertyFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF
)
:
    internalField_(iF),
    i_(-1)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::UPtrList<Foam::diameterModels::sizeGroup>&
Foam::secondaryPropertyFvScalarFieldSource::fis() const
{
    return
        refCast<const diameterModels::velocityGroup>
        (
            internalField_
           .db()
           .lookupObject<phaseSystem>(phaseSystem::propertiesName)
           .phases()[internalField_.group()]
           .diameter()
        ).popBal().sizeGroups();
}


Foam::label Foam::secondaryPropertyFvScalarFieldSource::i() const
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
            << "Could not determine the size-group index for "
            << "secondary property field " << internalField_.name()
            << exit(FatalError);
    }

    // Convert the trailing digits to an integer
    i_ = atoi(member(memberChari, member.size()).c_str());

    // Check that the index is in bounds
    if (i_ < 0 || i_ >= fis().size())
    {
        FatalErrorInFunction
            << "Size-group index for secondary property field "
            << internalField_.name() << " is out of range"
            << exit(FatalError);
    }

    return i_;
}


const Foam::diameterModels::sizeGroup&
Foam::secondaryPropertyFvScalarFieldSource::fi(const label deltai) const
{
    const UPtrList<diameterModels::sizeGroup>& fis = this->fis();

    const label i = this->i() + deltai;

    return
        i < 0 || i > fis.size() - 1
      ? NullObjectRef<diameterModels::sizeGroup>()
      : fis[i];
}


// ************************************************************************* //
