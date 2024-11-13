/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "rigidBodyPropellerDisk.H"
#include "motionSolver_fvMeshMover.H"
#include "motionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(rigidBodyPropellerDisk, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        rigidBodyPropellerDisk,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::rigidBodyPropellerDisk::readCoeffs(const dictionary& dict)
{
    propellerDisk::readCoeffs(dict);

    if (!dict.found("centre"))
    {
        FatalIOErrorInFunction(dict)
            << "Required entry 'centre' not found"
            << exit(FatalIOError);
    }

    dict.lookup("body") >> body_;
    bodyID_ = motion_.bodyIndex(body_);
}


Foam::vector Foam::fv::rigidBodyPropellerDisk::centre() const
{
    return motion_.p(bodyID_, centre_);
}


Foam::vector Foam::fv::rigidBodyPropellerDisk::normal() const
{
    return motion_.d(bodyID_, normal_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::rigidBodyPropellerDisk::rigidBodyPropellerDisk
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    propellerDisk(name, modelType, mesh, dict),
    motion_
    (
        refCast<const RBD::rigidBodyMotion>
        (
            refCast<const fvMeshMovers::motionSolver>(mesh.mover()).motion()
        )
    )
{
    readCoeffs(coeffs(dict));
}


// ************************************************************************* //
