/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

#include "functionalRigidBodyMeshMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(functionalRigidBodyMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        functionalRigidBodyMeshMotion,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::List<Foam::septernion>
Foam::functionalRigidBodyMeshMotion::transforms0() const
{
    if (SBMFs_.size() == 1)
    {
        return List<septernion>(1, SBMFs_[0].transformation());
    }
    else
    {
        List<septernion> transforms0(SBMFs_.size());

        forAll(SBMFs_, i)
        {
            transforms0[i] = SBMFs_[i].transformation();
        }

        return transforms0;
    }
}


void Foam::functionalRigidBodyMeshMotion::moveBodies()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionalRigidBodyMeshMotion::functionalRigidBodyMeshMotion
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    multiRigidBodyMeshMotion(name, mesh, dict)
{
    SBMFs_.setSize(bodyMeshes_.size());

    if (dict.isDict("bodies"))
    {
        const dictionary& bodiesDict = dict.subDict("bodies");

        label i = 0;
        forAllConstIter(IDLList<entry>, bodiesDict, iter)
        {
            const dictionary& bodyDict = iter().dict();

            SBMFs_.set
            (
                i++,
                solidBodyMotionFunction::New(bodyDict, mesh.time(), "function")
            );
        }
    }
    else
    {
        SBMFs_.set
        (
            0,
            solidBodyMotionFunction::New(dict, mesh.time(), "function")
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionalRigidBodyMeshMotion::~functionalRigidBodyMeshMotion()
{}


// ************************************************************************* //
