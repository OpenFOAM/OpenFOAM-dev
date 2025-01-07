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
#include "coupled.H"
#include "massive.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(cloud, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::cloud::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    eqn += coupledCloud_.carrierEqn(field);
}


template<class Type>
void Foam::fv::cloud::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    eqn += coupledCloud_.carrierEqn(field);
}


template<class Type>
void Foam::fv::cloud::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    eqn += coupledCloud_.carrierEqn(field);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::cloud::~cloud()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::fv::cloud::addsSupToField(const word& fieldName) const
{
    const word group = IOobject::group(fieldName);

    return
        cloud_.LagrangianModels().addsSupToField
        (
            IOobject::groupName(clouds::shaped::vName, group)
        )
     || cloud_.LagrangianModels().addsSupToField
        (
            IOobject::groupName(clouds::massive::mName, group)
        )
     || cloud_.LagrangianModels().addsSupToField(fieldName);
}


void Foam::fv::cloud::addSup(fvMatrix<scalar>& eqn) const
{
    NotImplemented;
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_FIELD_SUP, fv::cloud)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP, fv::cloud)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP, fv::cloud)


void Foam::fv::cloud::correct()
{
    cloudPtr_->solve();
}


void Foam::fv::cloud::preUpdateMesh()
{
    cloudPtr_->storePosition();
}


bool Foam::fv::cloud::movePoints()
{
    return true;
}


void Foam::fv::cloud::topoChange(const polyTopoChangeMap& map)
{
    cloudPtr_->topoChange(map);
}


void Foam::fv::cloud::mapMesh(const polyMeshMap& map)
{
    cloudPtr_->mapMesh(map);
}


void Foam::fv::cloud::distribute(const polyDistributionMap& map)
{
    cloudPtr_->distribute(map);
}


// ************************************************************************* //
