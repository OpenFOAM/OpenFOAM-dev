/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "Usher.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace packingDispersionModels
{
    defineTypeNameAndDebug(Usher, 0);
    addToRunTimeSelectionTable
    (
        packingDispersionModel,
        Usher,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T>
inline auto Foam::packingDispersionModels::Usher::sigma1(const T& alphad) const
{
    return sigma01_*pow
    (
        (alphad - alphaGel_)
       /(alphacp_ - alphad)*(b1_ + alphad - alphaGel_),
        n1_
    );
}


template<class T>
inline auto Foam::packingDispersionModels::Usher::sigmaPrime1
(
    const T& alphad
) const
{
    return n1_*sigma01_*pow
    (
        (alphad - alphaGel_)
       /(alphacp_ - alphad)*(b1_ + alphad - alphaGel_),
        n1_ + 1
    )
   *(b1_*(alphacp_ - alphaGel_)/sqr(alphad - alphaGel_) + 1);
}


template<class T>
inline auto Foam::packingDispersionModels::Usher::sigma2(const T& alphad) const
{
    return sigma02_*pow
    (
        (alphad - alphag_)
       /(alphaMax_ - alphad)*(b2_ + alphad - alphag_),
        n2_
    );
}


template<class T>
inline auto Foam::packingDispersionModels::Usher::sigmaPrime2
(
    const T& alphad
) const
{
    return n2_*sigma02_*pow
    (
        (alphad - alphag_)
       /(alphaMax_ - alphad)*(b2_ + alphad - alphag_),
        n2_ + 1
    )
   *(b2_*(alphaMax_ - alphag_)/sqr(alphad - alphag_) + 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingDispersionModels::Usher::Usher
(
    const dictionary& dict,
    const relativeVelocityModel& relativeVelocity
)
:
    packingDispersionModel(relativeVelocity),
    alphaGel_("alphaGel", dimless, dict),
    alphap_("alphap", dimless, dict),
    alphag_("alphap", dimless, dict),
    alphacp_("alphacp", dimless, dict),
    alphaMax_("alphaMax", dimless, dict),
    b1_("b1", dimless, dict),
    n1_("n1", dimless, dict),
    sigma01_("sigma01", sqr(dimVelocity), dict),
    b2_("b2", dimless, dict),
    n2_
    (
        "n2",
        sigmaPrime1(alphap_)/sigma1(alphap_)
       *(alphaMax_ - alphap_)*(alphap_ - alphag_)
       *(b2_ + alphap_ - alphag_)
       /((b2_*(alphaMax_ - alphag_) + sqr(alphap_ - alphag_)))
    ),
    sigma02_
    (
        "sigma02",
        sigma1(alphap_)
       *pow((alphaMax_ - alphap_)*(b2_ + alphap_ - alphag_), n2_)
       /pow(alphap_ - alphag_, n2_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingDispersionModels::Usher::~Usher()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::packingDispersionModels::Usher::sigmaPrime() const
{
    const volScalarField& alphad = mixture_.alphad();

    return
        pos(alphad - alphaGel_)*pos(alphap_ - alphad)*sigmaPrime1(alphad)
      + pos(alphad - alphap_)*pos(alphaMax_ - alphad)*sigmaPrime2(alphad);
}


// ************************************************************************* //
