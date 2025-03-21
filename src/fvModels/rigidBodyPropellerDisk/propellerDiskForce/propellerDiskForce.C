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

#include "propellerDiskForce.H"
#include "rigidBodyMeshMotion.H"
#include "fvModels.H"
#include "rigidBodyPropellerDisk.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(propellerDiskForce, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        propellerDiskForce,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::propellerDiskForce::propellerDiskForce
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

Foam::RBD::restraints::propellerDiskForce::~propellerDiskForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::propellerDiskForce::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{
    const rigidBodyMeshMotion& mover =
        refCast<const rigidBodyMeshMotion>(model_);

    const fvMesh& mesh = refCast<const fvMesh>(mover.mesh());

    // Lookup the fvModels for this mesh
    const fvModels& fvModels(fvModels::New(mesh));

    // Search for the rigidBodyPropellerDisk fvModel
    // for the body this restraint is attached to
    const fv::rigidBodyPropellerDisk* propPtr = nullptr;
    forAll(fvModels, i)
    {
        if (isType<fv::rigidBodyPropellerDisk>(fvModels[i]))
        {
            const fv::rigidBodyPropellerDisk& prop
            (
                refCast<const fv::rigidBodyPropellerDisk>(fvModels[i])
            );

            if (bodyIndex_ == prop.bodyID())
            {
                propPtr = &prop;
            }
        }
    }

    if (!propPtr)
    {
        FatalErrorInFunction
            << "Cannot find " << fv::rigidBodyPropellerDisk::typeName
            << " fvModel for body "
            << model_.name(bodyIndex_)
            << exit(FatalError);
    }

    // Lookup the propeller force and moment
    const vector force(propPtr->force());
    const vector moment
    (
        propPtr->moment()
      + ((propPtr->centre() - model_.X0(bodyIndex_).r()) ^ force)
    );

    if (model_.debug)
    {
        Info<< " location " << propPtr->centre()
            << " force " << force
            << " moment " << moment
            << endl;
    }

    // Accumulate the force for the restrained body
    fx[masterBodyIndex_] += spatialVector(moment, force);
}


bool Foam::RBD::restraints::propellerDiskForce::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    return true;
}


void Foam::RBD::restraints::propellerDiskForce::write
(
    Ostream& os
) const
{
    restraint::write(os);
}


// ************************************************************************* //
