/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "fvmSup.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::MomentumCloud<CloudType>&
Foam::MomentumCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::MomentumCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::MomentumCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::MomentumCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::MomentumCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::MomentumCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::MomentumCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::MomentumCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::MomentumCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::MomentumCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::MomentumCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::MomentumCloud<CloudType>::mu() const
{
    return mu_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::MomentumCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
inline Foam::scalar Foam::MomentumCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::MomentumCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}


template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::MomentumCloud<CloudType>::forceType&
Foam::MomentumCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::MomentumCloud<CloudType>::forceType&
Foam::MomentumCloud<CloudType>::forces()
{
    return forces_;
}


template<class CloudType>
inline typename Foam::MomentumCloud<CloudType>::functionType&
Foam::MomentumCloud<CloudType>::functions()
{
    return functions_;
}


template<class CloudType>
inline const Foam::InjectionModelList<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::injectors() const
{
    return injectors_;
}


template<class CloudType>
inline Foam::InjectionModelList<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::injectors()
{
    return injectors_;
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::dispersion() const
{
    return dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::dispersion()
{
    return dispersionModel_();
}


template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::patchInteraction() const
{
    return patchInteractionModel_;
}


template<class CloudType>
inline Foam::PatchInteractionModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::patchInteraction()
{
    return patchInteractionModel_();
}


template<class CloudType>
inline const Foam::StochasticCollisionModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::stochasticCollision() const
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline Foam::StochasticCollisionModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::stochasticCollision()
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::surfaceFilm() const
{
    return filmModel_();
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::MomentumCloud<CloudType>>&
Foam::MomentumCloud<CloudType>::surfaceFilm()
{
    return filmModel_();
}


template<class CloudType>
inline const Foam::integrationScheme&
Foam::MomentumCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class CloudType>
inline Foam::label Foam::MomentumCloud<CloudType>::nParcels() const
{
    return this->size();
}


template<class CloudType>
inline Foam::scalar Foam::MomentumCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::MomentumCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::MomentumCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::MomentumCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, vSmall);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::MomentumCloud<CloudType>::Dmin() const
{
    scalar d = vGreat;

    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        d = min(d, iter().d());
    }

    return returnReduce(d, minOp<scalar>());
}


template<class CloudType>
inline Foam::scalar Foam::MomentumCloud<CloudType>::Dmax() const
{
    scalar d = -vGreat;

    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        d = max(d, iter().d());
    }

    return returnReduce(d, maxOp<scalar>());
}


template<class CloudType>
inline Foam::Random& Foam::MomentumCloud<CloudType>::rndGen() const
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*>>&
Foam::MomentumCloud<CloudType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline const Foam::scalarField&
Foam::MomentumCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::MomentumCloud<CloudType>::UTrans() const
{
    return UTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::MomentumCloud<CloudType>::UTransRef()
{
    return UTrans_();
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::MomentumCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::MomentumCloud<CloudType>::UCoeffRef()
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::MomentumCloud<CloudType>::SU(const volVectorField& U) const
{
    if (debug)
    {
        Info<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            const volScalarField::Internal
                Vdt(this->mesh().V()*this->db().time().deltaT());

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm.ref();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::MomentumCloud<CloudType>::vDotSweep() const
{
    tmp<volScalarField> tvDotSweep
    (
        volScalarField::New
        (
            this->name() + ":vDotSweep",
            this->mesh(),
            dimensionedScalar(dimless/dimTime, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& vDotSweep = tvDotSweep.ref();
    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        vDotSweep[celli] += p.nParticle()*p.areaP()*mag(p.U() - U_[celli]);
    }

    vDotSweep.primitiveFieldRef() /= this->mesh().V();
    vDotSweep.correctBoundaryConditions();

    return tvDotSweep;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::MomentumCloud<CloudType>::alpha() const
{
    tmp<volScalarField> ttheta
    (
        volScalarField::New
        (
            this->name() + ":alpha",
            this->mesh(),
            dimensionedScalar(dimless, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& alpha = ttheta.ref();
    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        alpha[celli] += p.nParticle()*p.volume();
    }

    alpha.primitiveFieldRef() /= this->mesh().V();
    alpha.correctBoundaryConditions();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::MomentumCloud<CloudType>::Y() const
{
    tmp<volScalarField> tY
    (
        volScalarField::New
        (
            this->name() + ":Y",
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    scalarField& Y = tY.ref().primitiveFieldRef();
    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        Y[celli] += p.nParticle()*p.mass();
    }

    Y /= (this->mesh().V()*rho_);

    return tY;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::MomentumCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        volScalarField::New
        (
            this->name() + ":rhoEff",
            this->mesh(),
            dimensionedScalar(dimDensity, 0)
        )
    );

    scalarField& rhoEff = trhoEff.ref().primitiveFieldRef();
    forAllConstIter(typename MomentumCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        rhoEff[celli] += p.nParticle()*p.mass();
    }

    rhoEff /= this->mesh().V();

    return trhoEff;
}


// ************************************************************************* //
