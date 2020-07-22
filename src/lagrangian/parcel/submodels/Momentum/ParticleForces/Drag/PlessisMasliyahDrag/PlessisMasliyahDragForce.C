/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PlessisMasliyahDragForce<CloudType>::PlessisMasliyahDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    DenseDragForce<CloudType>(owner, mesh, dict, typeName)
{}


template<class CloudType>
Foam::PlessisMasliyahDragForce<CloudType>::PlessisMasliyahDragForce
(
    const PlessisMasliyahDragForce<CloudType>& df
)
:
    DenseDragForce<CloudType>(df)
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
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar alphac =
        this->alphacInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices()
        );

    const scalar cbrtAlphap = pow(1 - alphac, 1.0/3.0);

    const scalar A =
        26.8*pow3(alphac)
       /(
            sqr(cbrtAlphap)
           *(1 - cbrtAlphap)
           *sqr(1 - sqr(cbrtAlphap))
          + small
        );

    const scalar B =
        sqr(alphac)
       /sqr(1 - sqr(cbrtAlphap));

    return forceSuSp
    (
        Zero,
        (mass/p.rho())
       *(A*(1 - alphac)/alphac + B*Re)*muc/(alphac*sqr(p.d()))
    );
}


// ************************************************************************* //
