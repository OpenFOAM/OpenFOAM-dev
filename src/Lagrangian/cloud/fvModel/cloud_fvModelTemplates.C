/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cloud_fvModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::cloud::cloud
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const cloud::Cloud<Type>&
)
:
    fvModel(name, modelType, mesh, dict),
    cloudPtr_(new Type(mesh, name, Foam::cloud::contextType::fvModel)),
    cloud_(cloudPtr_()),
    coupledCloud_(refCast<const clouds::coupled>(cloud_))
{
    // Ensure LagrangianModels are constructed before time is incremented
    cloudPtr_->LagrangianModels();
}


template<class Type>
Foam::fv::Cloud<Type>::Cloud
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    cloud(name, modelType, mesh, dict, cloud::Cloud<Type>())
{}


// ************************************************************************* //
