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

#include "compositeJoint.H"
#include "rigidBodyModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(composite, 0);

    addToRunTimeSelectionTable
    (
        joint,
        composite,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::RBD::joints::composite::setLastJoint()
{
    last().joint::operator=(*this);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::composite::composite(const PtrList<joint>& joints)
:
    PtrList<joint>(joints),
    joint(last())
{}


Foam::RBD::joints::composite::composite
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    PtrList<joint>(dict.lookup("joints"), joint::iNew(model)),
    joint(last())
{}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::composite::clone() const
{
    return autoPtr<joint>(new composite(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::composite::~composite()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::composite::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    last().jcalc(J, state);
}


void Foam::RBD::joints::composite::write(Ostream& os) const
{
    joint::write(os);
    writeKeyword(os, "joints");
    os << static_cast<const PtrList<joint>&>(*this);
}


// ************************************************************************* //
