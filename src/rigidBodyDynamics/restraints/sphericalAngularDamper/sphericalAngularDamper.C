/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "sphericalAngularDamper.H"
#include "rigidBodyModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(sphericalAngularDamper, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        sphericalAngularDamper,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::sphericalAngularDamper::sphericalAngularDamper
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraints::sphericalAngularDamper::~sphericalAngularDamper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::sphericalAngularDamper::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{
    vector moment = -coeff_*model_.v(model_.master(bodyID_)).w();

    if (model_.debug)
    {
        Info<< " moment " << moment << endl;
    }

    // Accumulate the force for the restrained body
    fx[bodyIndex_] += model_.X0(bodyID_).T() & spatialVector(moment, Zero);
}


bool Foam::RBD::restraints::sphericalAngularDamper::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    coeffs_.lookup("coeff") >> coeff_;

    return true;
}


void Foam::RBD::restraints::sphericalAngularDamper::write
(
    Ostream& os
) const
{
    restraint::write(os);

    writeEntry(os, "coeff", coeff_);
}


// ************************************************************************* //
