/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "Stokes2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Stokes2, 0);
    addToRunTimeSelectionTable(waveModel, Stokes2, dictionary);
}
}


// * * * * * * * * * * Static Protected Member Functions  * * * * * * ** * * //

Foam::scalar Foam::waveModels::Stokes2::celerity(const AiryCoeffs& coeffs)
{
    static const scalar kdGreat = log(great);
    const scalar kd = min(max(coeffs.k()*coeffs.depth, - kdGreat), kdGreat);
    const scalar ka = coeffs.k()*coeffs.amplitude;

    const scalar S = coeffs.deep() ? 0 : 1/cosh(2*kd);

    const scalar C0 = coeffs.celerity();
    const scalar C2ByC0 = (2 + 7*sqr(S))/4/sqr(1 - S);

    if (debug)
    {
        Info<< "C2 = " << C2ByC0*C0 << endl;
    }

    return Airy::celerity(coeffs) + sqr(ka)*C2ByC0*C0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Stokes2::Stokes2
(
    const dictionary& dict,
    const scalar g,
    const word& modelName,
    scalar (*celerityPtr)(const AiryCoeffs&)
)
:
    Airy(dict, g, modelName, celerityPtr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::Stokes2::~Stokes2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Stokes2::celerity() const
{
    return celerity(coeffs());
}


Foam::tmp<Foam::scalarField> Foam::waveModels::Stokes2::elevation
(
    const scalar t,
    const scalarField& x
) const
{
    const AiryCoeffs coeffs = this->coeffs();

    static const scalar kdGreat = log(great);
    const scalar kd = min(max(coeffs.k()*depth(), - kdGreat), kdGreat);
    const scalar ka = coeffs.k()*amplitude(t);

    const scalar T = coeffs.deep() ? 1 : tanh(kd);

    const scalar B22 = (3/sqr(T) - 1)/T/4;

    if (debug)
    {
        Info<< "B22 = " << B22 << endl;
    }

    return
        Airy::elevation(t, x)
      + (1/coeffs.k())*sqr(ka)*B22*cos(2*coeffs.angle(phase(), t, x));
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Stokes2::velocity
(
    const scalar t,
    const vector2DField& xz
) const
{
    const AiryCoeffs coeffs = this->coeffs();

    static const scalar kdGreat = log(great);
    const scalar kd = min(max(coeffs.k()*depth(), - kdGreat), kdGreat);
    const scalar ka = coeffs.k()*amplitude(t);

    const scalar A22ByA11 = coeffs.deep() ? 0 : 0.375/pow3(sinh(kd));

    if (debug)
    {
        const scalar A11 = 1/sinh(kd);
        Info<< "A22 = " << A22ByA11*A11 << endl;
    }

    return
        Airy::velocity(t, xz)
      + Airy::celerity()*sqr(ka)*A22ByA11*coeffs.vi(2, phase(), t, xz);
}


// ************************************************************************* //
