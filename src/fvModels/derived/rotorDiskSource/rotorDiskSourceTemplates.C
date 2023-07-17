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

#include "rotorDiskSource.H"
#include "volFields.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::rotorDiskSource::calculate
(
    const RhoFieldType& rho,
    const vectorField& U,
    const scalarField& thetag,
    vectorField& force,
    const bool divideVolume,
    const bool output
) const
{
    const scalarField& V = mesh().V();

    // Logging info
    scalar dragEff = 0;
    scalar liftEff = 0;
    scalar AOAmin = great;
    scalar AOAmax = -great;
    scalar powerEff = 0;

    const labelUList cells = set_.cells();

    forAll(cells, i)
    {
        if (area_[i] > rootVSmall)
        {
            const label celli = cells[i];

            const scalar radius = x_[i].x();

            // Transform velocity into local cylindrical reference frame
            vector Uc = cylindrical_->invTransform(U[celli], i);
            // Uc.x(): radial direction.
            // Uc.y(): drag direction.
            // Uc.z(): lift / thrust direction.

            // Transform velocity into local coning system
            Uc = R_[i] & Uc;

            // Set radial component of velocity to zero
            Uc.x() = 0;

            // Set blade normal component of velocity
            Uc.y() = radius*omega_ - Uc.y();

            // Determine blade data for this radius
            // i2 = index of upper radius bound data point in blade list
            scalar twist = 0;
            scalar chord = 0;
            label i1 = -1;
            label i2 = -1;
            scalar invDr = 0;
            blade_.interpolate(radius, twist, chord, i1, i2, invDr);

            const scalar alphaGeom = thetag[i] + twist;

            // Effective angle of attack
            const int rotationSign = sign(omega_);
            const scalar alphaEff =
                alphaGeom - atan2(-Uc.z(), rotationSign*Uc.y());

            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            // Determine profile data for this radius and angle of attack
            const label profile1 = blade_.profileID()[i1];
            const label profile2 = blade_.profileID()[i2];

            scalar Cd1 = 0;
            scalar Cl1 = 0;
            profiles_[profile1].Cdl(alphaEff, Cd1, Cl1);

            scalar Cd2 = 0;
            scalar Cl2 = 0;
            profiles_[profile2].Cdl(alphaEff, Cd2, Cl2);

            const scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
            const scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

            // Apply tip effect for blade lift
            const scalar tipFactor = neg(radius/rMax_ - tipEffect_);

            // Calculate forces perpendicular to blade
            const scalar pDyn = 0.5*rho[celli]*magSqr(Uc);

            const scalar f =
                pDyn*chord*nBlades_*area_[i]/radius/mathematical::twoPi;

            vector localForce = vector(0, rotationSign*-f*Cd, tipFactor*f*Cl);

            // Accumulate forces
            dragEff += rhoRef_*localForce.y();
            liftEff += rhoRef_*localForce.z();
            powerEff += rhoRef_*localForce.y()*radius*omega_;

            // Transform force from local coning system into rotor cylindrical
            localForce = invR_[i] & localForce;

            // Transform force into global Cartesian co-ordinate system
            force[celli] = cylindrical_->transform(localForce, i);

            if (divideVolume)
            {
                force[celli] /= V[celli];
            }
        }
    }

    if (output)
    {
        reduce(AOAmin, minOp<scalar>());
        reduce(AOAmax, maxOp<scalar>());
        reduce(dragEff, sumOp<scalar>());
        reduce(liftEff, sumOp<scalar>());

        Info<< type() << " output:" << nl
            << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
            << radToDeg(AOAmax) << nl
            << "    Effective power = " << powerEff << nl
            << "    Effective drag = " << dragEff << nl
            << "    Effective lift = " << liftEff << endl;
    }
}


template<class Type>
void Foam::fv::rotorDiskSource::writeField
(
    const word& name,
    const List<Type>& values
) const
{
    tmp<VolInternalField<Type>> tfield
    (
        new VolField<Type>
        (
            IOobject
            (
                name,
                mesh().time().name(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensioned<Type>("zero", dimless, Zero)
        )
    );

    Field<Type>& field = tfield.ref();

    const labelUList cells = set_.cells();

    if (cells.size() != values.size())
    {
        FatalErrorInFunction << abort(FatalError);
    }

    forAll(cells, i)
    {
        field[cells[i]] = values[i];
    }

    tfield().write();
}


// ************************************************************************* //
