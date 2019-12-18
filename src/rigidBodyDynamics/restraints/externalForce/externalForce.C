/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "externalForce.H"
#include "rigidBodyModel.H"
#include "rigidBodyModelState.H"
#include "OneConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(externalForce, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        externalForce,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::externalForce::externalForce
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model),
    externalForce_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraints::externalForce::~externalForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::externalForce::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{
    const vector force = externalForce_().value(state.t());
    const vector moment(location_ ^ force);

    if (model_.debug)
    {
        Info<< " location " << location_
            << " force " << force
            << " moment " << moment
            << endl;
    }

    // Accumulate the force for the restrained body
    fx[bodyIndex_] += spatialVector(moment, force);
}


bool Foam::RBD::restraints::externalForce::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    coeffs_.lookup("location") >> location_;

    externalForce_ = Function1<vector>::New("force", coeffs_);

    return true;
}


void Foam::RBD::restraints::externalForce::write
(
    Ostream& os
) const
{
    restraint::write(os);

    writeEntry(os, "location", location_);

    writeEntry(os, externalForce_());
}


// ************************************************************************* //
