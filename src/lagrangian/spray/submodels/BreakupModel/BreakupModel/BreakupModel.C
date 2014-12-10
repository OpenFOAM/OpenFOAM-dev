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

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner),
    solveOscillationEq_(false),
    y0_(0.0),
    yDot0_(0.0),
    TABComega_(0.0),
    TABCmu_(0.0),
    TABWeCrit_(0.0)
{}


template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    const BreakupModel<CloudType>& bum
)
:
    CloudSubModelBase<CloudType>(bum),
    solveOscillationEq_(bum.solveOscillationEq_),
    y0_(bum.y0_),
    yDot0_(bum.yDot0_),
    TABComega_(bum.TABComega_),
    TABCmu_(bum.TABCmu_),
    TABWeCrit_(bum.TABWeCrit_)
{}


template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    bool solveOscillationEq
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    solveOscillationEq_(solveOscillationEq),
    y0_(0.0),
    yDot0_(0.0),
    TABComega_(0.0),
    TABCmu_(0.0),
    TABWeCrit_(0.0)
{
    if (solveOscillationEq_)
    {
        const dictionary coeffs(dict.subDict("TABCoeffs"));
        y0_ = coeffs.template lookupOrDefault<scalar>("y0", 0.0);
        yDot0_ = coeffs.template lookupOrDefault<scalar>("yDot0", 0.0);
        TABComega_ = coeffs.template lookupOrDefault<scalar>("Comega", 8.0);
        TABCmu_ = coeffs.template lookupOrDefault<scalar>("Cmu", 10.0);
        TABWeCrit_ = coeffs.template lookupOrDefault<scalar>("WeCrit", 12.0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BreakupModel<CloudType>::~BreakupModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::BreakupModel<CloudType>::update
(
    const scalar dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar d0,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const vector& U,
    const scalar rhoc,
    const scalar muc,
    const vector& Urel,
    const scalar Urmag,
    const scalar tMom,
    scalar& dChild,
    scalar& massChild
)
{
    notImplemented
    (
        "bool Foam::BreakupModel<CloudType>::update"
        "("
            "const scalar, "
            "const vector&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const vector&, "
            "const scalar, "
            "const scalar, "
            "const vector&, "
            "const scalar, "
            "const scalar, "
            "scalar&, "
            "scalar&"
        ");"
    );

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BreakupModelNew.C"

// ************************************************************************* //

