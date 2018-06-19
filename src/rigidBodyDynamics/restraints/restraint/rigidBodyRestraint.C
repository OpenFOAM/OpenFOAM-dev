/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "rigidBodyRestraint.H"
#include "rigidBodyModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(restraint, 0);
    defineRunTimeSelectionTable(restraint, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraint::restraint
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    name_(name),
    bodyID_(model.bodyID(dict.lookup("body"))),
    bodyIndex_(model.master(bodyID_)),
    coeffs_(dict),
    model_(model)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraint::~restraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::RBD::restraint::coeffDict() const
{
    return coeffs_;
}


bool Foam::RBD::restraint::read(const dictionary& dict)
{
    coeffs_ = dict;
    return true;
}


void Foam::RBD::restraint::write(Ostream& os) const
{
    os.writeKeyword("type")
        << type() << token::END_STATEMENT << nl;
    os.writeKeyword("body")
        << model_.name(bodyID_) << token::END_STATEMENT << nl;
}


// ************************************************************************* //
