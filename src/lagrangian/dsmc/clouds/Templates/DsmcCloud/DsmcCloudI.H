/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::word& Foam::DsmcCloud<ParcelType>::cloudName() const
{
    return cloudName_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::DsmcCloud<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::DsmcCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::List<Foam::word>&
Foam::DsmcCloud<ParcelType>::typeIdList() const
{
    return typeIdList_;
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline const Foam::List<Foam::DynamicList<ParcelType*> >&
Foam::DsmcCloud<ParcelType>::cellOccupancy() const
{
    return cellOccupancy_;
}


template<class ParcelType>
inline Foam::volScalarField& Foam::DsmcCloud<ParcelType>::sigmaTcRMax()
{
    return sigmaTcRMax_;
}


template<class ParcelType>
inline Foam::scalarField&
Foam::DsmcCloud<ParcelType>::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


template<class ParcelType>
inline const Foam::List<typename ParcelType::constantProperties>&
Foam::DsmcCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::DsmcCloud<ParcelType>::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorIn("Foam::DsmcCloud<ParcelType>::constProps(label typeId)")
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}


template<class ParcelType>
inline Foam::Random& Foam::DsmcCloud<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline Foam::volScalarField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::qBF()
{
    return q_.boundaryField();
}


template<class ParcelType>
inline Foam::volVectorField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::fDBF()
{
    return fD_.boundaryField();
}


template<class ParcelType>
inline Foam::volScalarField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::rhoNBF()
{
    return rhoN_.boundaryField();
}


template<class ParcelType>
inline Foam::volScalarField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::rhoMBF()
{
    return rhoM_.boundaryField();
}


template<class ParcelType>
inline Foam::volScalarField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::linearKEBF()
{
    return linearKE_.boundaryField();
}


template<class ParcelType>
inline Foam::volScalarField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::internalEBF()
{
    return internalE_.boundaryField();
}


template<class ParcelType>
inline Foam::volScalarField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::iDofBF()
{
    return iDof_.boundaryField();
}


template<class ParcelType>
inline Foam::volVectorField::GeometricBoundaryField&
Foam::DsmcCloud<ParcelType>::momentumBF()
{
    return momentum_.boundaryField();
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DsmcCloud<ParcelType>::boundaryT() const
{
    return boundaryT_;
}


template<class ParcelType>
inline const Foam::volVectorField&
Foam::DsmcCloud<ParcelType>::boundaryU() const
{
    return boundaryU_;
}


template<class ParcelType>
inline const Foam::BinaryCollisionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::binaryCollision() const
{
    return binaryCollisionModel_;
}


template<class ParcelType>
inline Foam::BinaryCollisionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::binaryCollision()
{
    return binaryCollisionModel_();
}


template<class ParcelType>
inline const Foam::WallInteractionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::wallInteraction() const
{
    return wallInteractionModel_;
}


template<class ParcelType>
inline Foam::WallInteractionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::wallInteraction()
{
    return wallInteractionModel_();
}


template<class ParcelType>
inline const Foam::InflowBoundaryModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::inflowBoundary() const
{
    return inflowBoundaryModel_;
}


template<class ParcelType>
inline Foam::InflowBoundaryModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::inflowBoundary()
{
    return inflowBoundaryModel_();
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysMass += cP.mass();
    }

    return nParticle_*sysMass;
}


template<class ParcelType>
inline Foam::vector Foam::DsmcCloud<ParcelType>::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum += cP.mass()*p.U();
    }

    return nParticle_*linearMomentum;
}


template<class ParcelType>
inline Foam::scalar
Foam::DsmcCloud<ParcelType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U());
    }

    return nParticle_*linearKineticEnergy;
}


template<class ParcelType>
inline Foam::scalar
Foam::DsmcCloud<ParcelType>::internalEnergyOfSystem() const
{
    scalar internalEnergy = 0.0;

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        internalEnergy += p.Ei();
    }

    return nParticle_*internalEnergy;
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


template<class ParcelType>
inline Foam::scalarField Foam::DsmcCloud<ParcelType>::maxwellianAverageSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
    return tfld();
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField Foam::DsmcCloud<ParcelType>::maxwellianRMSSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(3.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}


template<class ParcelType>
inline Foam::scalar
Foam::DsmcCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField
Foam::DsmcCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(2.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::DsmcCloud<ParcelType>::q() const
{
    return q_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DsmcCloud<ParcelType>::fD() const
{
    return fD_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DsmcCloud<ParcelType>::rhoN() const
{
    return rhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::DsmcCloud<ParcelType>::rhoM() const
{
    return rhoM_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DsmcCloud<ParcelType>::dsmcRhoN() const
{
    return dsmcRhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DsmcCloud<ParcelType>::linearKE() const
{
    return linearKE_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DsmcCloud<ParcelType>::internalE() const
{
    return internalE_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DsmcCloud<ParcelType>::iDof() const
{
    return iDof_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DsmcCloud<ParcelType>::momentum() const
{
    return momentum_;
}


template<class ParcelType>
inline void Foam::DsmcCloud<ParcelType>::clear()
{
    return IDLList<ParcelType>::clear();
}


// ************************************************************************* //
