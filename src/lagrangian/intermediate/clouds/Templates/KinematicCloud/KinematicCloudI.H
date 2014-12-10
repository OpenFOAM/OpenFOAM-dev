/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::KinematicCloud<CloudType>&
Foam::KinematicCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::KinematicCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::KinematicCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::KinematicCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::KinematicCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::KinematicCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::KinematicCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::KinematicCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicCloud<CloudType>::mu() const
{
    return mu_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::KinematicCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::KinematicCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}


template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::KinematicCloud<CloudType>::forceType&
Foam::KinematicCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::KinematicCloud<CloudType>::forceType&
Foam::KinematicCloud<CloudType>::forces()
{
    return forces_;
}


template<class CloudType>
inline typename Foam::KinematicCloud<CloudType>::functionType&
Foam::KinematicCloud<CloudType>::functions()
{
    return functions_;
}


template<class CloudType>
inline const Foam::InjectionModelList<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::injectors() const
{
    return injectors_;
}


template<class CloudType>
inline Foam::InjectionModelList<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::injectors()
{
    return injectors_;
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::dispersion() const
{
    return dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::dispersion()
{
    return dispersionModel_();
}


template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::patchInteraction() const
{
    return patchInteractionModel_;
}


template<class CloudType>
inline Foam::PatchInteractionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::patchInteraction()
{
    return patchInteractionModel_();
}


template<class CloudType>
inline const Foam::StochasticCollisionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::stochasticCollision() const
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline Foam::StochasticCollisionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::stochasticCollision()
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::surfaceFilm() const
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::surfaceFilm()
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline const Foam::vectorIntegrationScheme&
Foam::KinematicCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class CloudType>
inline Foam::label Foam::KinematicCloud<CloudType>::nParcels() const
{
    return this->size();
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::KinematicCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::Dmax() const
{
    scalar d = -GREAT;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::penetration
(
    const scalar fraction
) const
{
    if ((fraction < 0) || (fraction > 1))
    {
        FatalErrorIn
        (
            "inline Foam::scalar Foam::KinematicCloud<CloudType>::penetration"
            "("
                "const scalar"
            ") const"
        )
            << "fraction should be in the range 0 < fraction < 1"
            << exit(FatalError);
    }

    scalar distance = 0.0;

    const label nParcel = this->size();
    globalIndex globalParcels(nParcel);
    const label nParcelSum = globalParcels.size();

    if (nParcelSum == 0)
    {
        return distance;
    }

    // lists of parcels mass and distance from initial injection point
    List<List<scalar> > procMass(Pstream::nProcs());
    List<List<scalar> > procDist(Pstream::nProcs());

    List<scalar>& mass = procMass[Pstream::myProcNo()];
    List<scalar>& dist = procDist[Pstream::myProcNo()];

    mass.setSize(nParcel);
    dist.setSize(nParcel);

    label i = 0;
    scalar mSum = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        scalar m = p.nParticle()*p.mass();
        scalar d = mag(p.position() - p.position0());
        mSum += m;

        mass[i] = m;
        dist[i] = d;

        i++;
    }

    // calculate total mass across all processors
    reduce(mSum, sumOp<scalar>());
    Pstream::gatherList(procMass);
    Pstream::gatherList(procDist);

    if (Pstream::master())
    {
        // flatten the mass lists
        List<scalar> allMass(nParcelSum, 0.0);
        SortableList<scalar> allDist(nParcelSum, 0.0);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            SubList<scalar>
            (
                allMass,
                globalParcels.localSize(procI),
                globalParcels.offset(procI)
            ).assign(procMass[procI]);

            // flatten the distance list
            SubList<scalar>
            (
                allDist,
                globalParcels.localSize(procI),
                globalParcels.offset(procI)
            ).assign(procDist[procI]);
        }

        // sort allDist distances into ascending order
        // note: allMass masses are left unsorted
        allDist.sort();

        if (nParcelSum > 1)
        {
            const scalar mLimit = fraction*mSum;
            const labelList& indices = allDist.indices();

            if (mLimit > (mSum - allMass[indices.last()]))
            {
                distance = allDist.last();
            }
            else
            {
                // assuming that 'fraction' is generally closer to 1 than 0,
                // loop through in reverse distance order
                const scalar mThreshold = (1.0 - fraction)*mSum;
                scalar mCurrent = 0.0;
                label i0 = 0;

                forAllReverse(indices, i)
                {
                    label indI = indices[i];

                    mCurrent += allMass[indI];

                    if (mCurrent > mThreshold)
                    {
                        i0 = i;
                        break;
                    }
                }

                if (i0 == indices.size() - 1)
                {
                    distance = allDist.last();
                }
                else
                {
                    // linearly interpolate to determine distance
                    scalar alpha = (mCurrent - mThreshold)/allMass[indices[i0]];
                    distance =
                        allDist[i0] + alpha*(allDist[i0+1] - allDist[i0]);
                }
            }
        }
        else
        {
            distance = allDist.first();
        }
    }

    Pstream::scatter(distance);

    return distance;
}


template<class CloudType>
inline Foam::cachedRandom& Foam::KinematicCloud<CloudType>::rndGen()
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*> >&
Foam::KinematicCloud<CloudType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline const Foam::scalarField&
Foam::KinematicCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTrans()
{
    return UTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTrans() const
{
    return UTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UCoeff()
{
    return UCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::KinematicCloud<CloudType>::SU(volVectorField& U) const
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
            const DimensionedField<scalar, volMesh>
                Vdt(mesh_.V()*this->db().time().deltaT());

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::vDotSweep() const
{
    tmp<volScalarField> tvDotSweep
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":vDotSweep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& vDotSweep = tvDotSweep();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        vDotSweep[cellI] += p.nParticle()*p.areaP()*mag(p.U() - U_[cellI]);
    }

    vDotSweep.internalField() /= mesh_.V();
    vDotSweep.correctBoundaryConditions();

    return tvDotSweep;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& theta = ttheta();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        theta[cellI] += p.nParticle()*p.volume();
    }

    theta.internalField() /= mesh_.V();
    theta.correctBoundaryConditions();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        alpha[cellI] += p.nParticle()*p.mass();
    }

    alpha /= (mesh_.V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":rhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        rhoEff[cellI] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh_.V();

    return trhoEff;
}


// ************************************************************************* //
