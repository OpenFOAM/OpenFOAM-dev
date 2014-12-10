/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "PlessisMasliyahDragForce.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::PlessisMasliyahDragForce<CloudType>::CdRe
(
    const scalar Re
) const
{
    if (Re > 1000.0)
    {
        return 0.44*Re;
    }
    else
    {
        return 24.0*(1.0 + 0.15*pow(Re, 0.687));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PlessisMasliyahDragForce<CloudType>::PlessisMasliyahDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    alphac_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            this->coeffs().lookup("alphac")
        )
    )
{}


template<class CloudType>
Foam::PlessisMasliyahDragForce<CloudType>::PlessisMasliyahDragForce
(
    const PlessisMasliyahDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    alphac_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            this->coeffs().lookup("alphac")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PlessisMasliyahDragForce<CloudType>::~PlessisMasliyahDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::PlessisMasliyahDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    scalar alphac(alphac_[p.cell()]);

    scalar cbrtAlphap(pow(1.0 - alphac, 1.0/3.0));

    scalar A =
        26.8*pow3(alphac)
       /(
            sqr(cbrtAlphap)
           *(1.0 - cbrtAlphap)
           *sqr(1.0 - sqr(cbrtAlphap))
          + SMALL
        );

    scalar B =
        sqr(alphac)
       /sqr(1.0 - sqr(cbrtAlphap));

    return forceSuSp
    (
        vector::zero,
        (mass/p.rho())
       *(A*(1.0 - alphac)/alphac + B*Re)*muc/(alphac*sqr(p.d()))
    );
}


// ************************************************************************* //
