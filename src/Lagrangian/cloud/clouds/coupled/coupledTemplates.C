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

#include "coupled.H"
#include "CarrierEqn.H"
#include "LagrangianmDdt.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::CarrierEqn<Type>& Foam::clouds::coupled::carrierEqn
(
    const CloudDerivedField<Type>& psic,
    const LagrangianSubMesh& subMesh
)
{
    const word key = carrierNameToName(subMesh.complete(psic(subMesh).name()));

    typename HashPtrTable<CarrierEqn<Type>>::iterator iter =
        carrierEqns<Type>().find(key);

    if (iter != carrierEqns<Type>().end()) return **iter;

    CarrierEqn<Type>* ptr = new CarrierEqn<Type>(key, Uc.psi().mesh());
    carrierEqns<Type>().insert(key, ptr);

    return *ptr;
}


template<class Type>
Foam::CarrierEqn<Type>& Foam::clouds::coupled::carrierEqn
(
    const CarrierField<Type>& psic,
    const LagrangianSubMesh&
)
{
    typename HashPtrTable<CarrierEqn<Type>>::iterator iter =
        carrierEqns<Type>().find(psic.psi().name());

    if (iter != carrierEqns<Type>().end()) return **iter;

    CarrierEqn<Type>* ptr = new CarrierEqn<Type>(psic.psi());
    carrierEqns<Type>().insert(psic.psi().name(), ptr);

    return *ptr;
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class Fieldc>
bool Foam::clouds::coupled::reCalculateModified
(
    const LagrangianSubScalarSubField& vOrM,
    const Fieldc& fieldc
)
{
    if (cloud_.context == cloud::contextType::functionObject) return false;

    const LagrangianSubMesh& subMesh = vOrM.mesh();

    return Lagrangianm::initDdt(vOrM.dimensions(), fieldc(subMesh));
}


template<class Fieldc, class Scale>
void Foam::clouds::coupled::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final,
    const LagrangianSubScalarSubField& vOrM,
    const Fieldc& fieldc,
    const Scale& scale
)
{
    if (cloud_.context == cloud::contextType::functionObject || !final) return;

    const LagrangianModels& models = cloud_.LagrangianModels();
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianEqn<scalar> fieldcEqn
    (
        Lagrangianm::noDdt(deltaT, vOrM.dimensions(), fieldc(subMesh))
     ==
        models.sourceProxy(deltaT, vOrM, fieldc(subMesh))
    );

    carrierEqn(fieldc, subMesh) += scale*fieldcEqn;
}


template<class Type, class Fieldc, class Scale>
void Foam::clouds::coupled::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubSubField<Type>& field,
    const Fieldc& fieldc,
    const Scale& scale
)
{
    if (cloud_.context == cloud::contextType::functionObject || !final) return;

    const LagrangianModels& models = cloud_.LagrangianModels();
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianEqn<Type> fieldcEqn
    (
        Lagrangianm::noDdt(deltaT, vOrM.dimensions(), fieldc(subMesh))
     ==
        models.sourceProxy(deltaT, vOrM, field, fieldc(subMesh))
    );

    carrierEqn(fieldc, subMesh) += scale*fieldcEqn;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::clouds::coupled::hasCarrierEqn(const word& key) const
{
    return carrierEqns<Type>().found(key);
}


template<class Type>
const Foam::CarrierEqn<Type>&
Foam::clouds::coupled::carrierEqn(const word& key) const
{
    return *carrierEqns<Type>()[key];
}


template<class Type>
const Foam::HashTable<const Foam::CarrierEqn<Type>*>
Foam::clouds::coupled::carrierEqns(const word& key) const
{
    const word member = IOobject::member(key);

    HashTable<const CarrierEqn<Type>*> result;

    forAllConstIter
    (
        typename HashPtrTable<CarrierEqn<Type>>,
        carrierEqns<Type>(),
        iter
    )
    {
        if (IOobject::member(iter.key()) == member)
        {
            result.set(IOobject::group(iter.key()), iter());
        }
    }

    return result;
}


template<class Type>
bool Foam::clouds::coupled::hasCarrierEqn(const VolField<Type>& psi) const
{
    return hasCarrierEqn<Type>(psi.name());
}


template<class Type>
const Foam::CarrierEqn<Type>&
Foam::clouds::coupled::carrierEqn(const VolField<Type>& psi) const
{
    return carrierEqn<Type>(psi.name());
}


template<class Type>
const Foam::HashTable<const Foam::CarrierEqn<Type>*>
Foam::clouds::coupled::carrierEqns(const VolField<Type>& psi) const
{
    return carrierEqns<Type>(psi.name());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
