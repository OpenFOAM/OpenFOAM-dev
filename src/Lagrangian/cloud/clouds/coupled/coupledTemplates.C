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
#include "CarrierEqn.H"
#include "LagrangianmDdt.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class Type>
Foam::CarrierEqn<Type>& Foam::clouds::coupled::carrierEqn
(
    const CloudDerivedField<Type>& psic
)
{
    const word key = carrierNameToName(psic.name());

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
    const CarrierField<Type>& psic
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

template<class Type>
bool Foam::clouds::coupled::initPsicDdt
(
    const LagrangianSubScalarSubField& vOrM,
    const CloudDerivedField<Type>& psic
) const
{
    const LagrangianSubMesh& subMesh = vOrM.mesh();

    return Lagrangianm::initDdt(vOrM.dimensions(), psic(subMesh));
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::clouds::coupled::psicEqn
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubSubField<Type>& psi,
    const CloudDerivedField<Type>& psic
) const
{
    const LagrangianModels& models = cloud_.LagrangianModels();
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    return
        Lagrangianm::noDdt(deltaT, vOrM.dimensions(), psic(subMesh))
     ==
        models.sourceProxy(deltaT, vOrM, psi, psic(subMesh));
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
