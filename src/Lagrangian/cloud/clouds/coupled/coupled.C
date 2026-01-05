/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "coupled.H"
#include "physicalProperties.H"
#include "viscosity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupled, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::clouds::coupled::getNucVf() const
{
    const word nucName = IOobject::groupName("nu", carrierPhaseName());

    if (cloud_.mesh().mesh().foundObject<volScalarField>(nucName))
    {
        return cloud_.mesh().mesh().lookupObject<volScalarField>(nucName);
    }

    const word viscosityName =
        IOobject::groupName(physicalProperties::typeName, carrierPhaseName());

    if (cloud_.mesh().mesh().foundObject<viscosity>(viscosityName))
    {
        return cloud_.mesh().mesh().lookupObject<viscosity>(viscosityName).nu();
    }

    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::LagrangianSubScalarField> Foam::clouds::coupled::calcNuc
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    FatalErrorInFunction
        << "Could not determine the carrier viscosity"
        << exit(FatalError);

    return tmp<LagrangianSubScalarField>(nullptr);
}


#define ACCESS_CARRIER_EQNS(Type, nullArg)                                     \
namespace Foam                                                                 \
{                                                                              \
    namespace clouds                                                           \
    {                                                                          \
        template<>                                                             \
        HashPtrTable<CarrierEqn<Type>>& coupled::carrierEqns() const           \
        {                                                                      \
            return CAT3(carrier, CAPITALIZE(Type), Eqns_);                     \
        }                                                                      \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_CARRIER_EQNS)
#undef ACCESS_CARRIER_EQNS


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupled::updateCarrier()
{
    if (tnucVf_.isTmp())
    {
        tnucVf_.ref() = getNucVf();
    }
}


void Foam::clouds::coupled::clearCarrierEqns()
{
    #define CLEAR_TYPE_CARRIER_EQNS(Type, nullArg)                             \
        forAllIter                                                             \
        (                                                                      \
            HashPtrTable<CarrierEqn<Type>>,                                    \
            carrierEqns<Type>(),                                               \
            iter                                                               \
        )                                                                      \
        {                                                                      \
            iter()->clear();                                                   \
        }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_CARRIER_EQNS);
    #undef CLEAR_TYPE_CARRIER_EQNS
}


Foam::tmp<Foam::LagrangianEqn<Foam::scalar>> Foam::clouds::coupled::psicEqn
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const CloudDerivedField<scalar>& oneOrRhoc
) const
{
    const LagrangianModels& models = cloud_.LagrangianModels();
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    return
        Lagrangianm::noDdt(deltaT, vOrM.dimensions(), oneOrRhoc(subMesh))
     ==
        models.sourceProxy(deltaT, vOrM, oneOrRhoc(subMesh));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupled::coupled(const cloud& c, const dictionary& dict)
:
    carried(c, dict),
    cloud_(c),
    tnucVf_(getNucVf()),
    nuc
    (
        tnucVf_.valid()
      ? carrierField<scalar>(tnucVf_())
      : c.derivedField<scalar>(*this, &coupled::calcNuc)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupled::~coupled()
{}


// ************************************************************************* //
