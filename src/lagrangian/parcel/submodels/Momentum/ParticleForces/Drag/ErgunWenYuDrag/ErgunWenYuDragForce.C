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

#include "ErgunWenYuDragForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ErgunWenYuDragForce<CloudType>::ErgunWenYuDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    WenYuDragForce<CloudType>(owner, mesh, dict, typeName)
{}


template<class CloudType>
Foam::ErgunWenYuDragForce<CloudType>::ErgunWenYuDragForce
(
    const ErgunWenYuDragForce<CloudType>& df
)
:
    WenYuDragForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ErgunWenYuDragForce<CloudType>::~ErgunWenYuDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::ErgunWenYuDragForce<CloudType>::calcCoupled
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

    if (alphac < 0.8)
    {
        return forceSuSp
        (
            Zero,
            mass/p.rho()*(150*(1 - alphac)/alphac + 1.75*Re)
           *muc/(alphac*sqr(p.d()))
        );
    }
    else
    {
        return WenYuDragForce<CloudType>::calcCoupled(p, td, dt, mass, Re, muc);
    }
}


// ************************************************************************* //
