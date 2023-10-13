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

#include "AiryCoeffs.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

Foam::scalar Foam::waveModels::AiryCoeffs::calcLength
(
    const scalar depth,
    const scalar amplitude,
    const scalar period,
    const scalar g,
    scalar (*celerityPtr)(const AiryCoeffs&)
)
{
    scalar length0 = 0;
    scalar length1 = 1.5*g*sqr(period)/(2*constant::mathematical::pi);

    // Bisect down to round off error to find the wave length
    for (label i = 0; i < ceil(std::log2(1/small)); ++ i)
    {
        const scalar length = (length0 + length1)/2;

        const scalar t =
            length/celerityPtr(AiryCoeffs(depth, amplitude, length, g));

        (t < period ? length0 : length1) = length;
    }

    return (length0 + length1)/2;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::AiryCoeffs::AiryCoeffs
(
    const scalar depth,
    const scalar amplitude,
    const scalar length,
    const scalar g
)
:
    depth(depth),
    amplitude(amplitude),
    length(length),
    g(g)
{}


Foam::waveModels::AiryCoeffs::AiryCoeffs
(
    const scalar depth,
    const scalar amplitude,
    const scalar period,
    const scalar g,
    scalar (*celerityPtr)(const AiryCoeffs&)
)
:
    depth(depth),
    amplitude(amplitude),
    length(calcLength(depth, amplitude, period, g, celerityPtr)),
    g(g)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::AiryCoeffs::k() const
{
    return 2*Foam::constant::mathematical::pi/length;
}


bool Foam::waveModels::AiryCoeffs::deep() const
{
    return depth*k() > log(great);
}


Foam::scalar Foam::waveModels::AiryCoeffs::celerity(const AiryCoeffs& coeffs)
{
    return sqrt(coeffs.g/coeffs.k()*tanh(coeffs.k()*coeffs.depth));
}


Foam::scalar Foam::waveModels::AiryCoeffs::celerity() const
{
    return celerity(*this);
}


Foam::tmp<Foam::scalarField> Foam::waveModels::AiryCoeffs::angle
(
    const scalar phase,
    const scalar t,
    const scalarField& x
) const
{
    return phase + k()*(x - celerity()*t);
}


Foam::tmp<Foam::scalarField> Foam::waveModels::AiryCoeffs::elevation
(
    const scalar phase,
    const scalar t,
    const scalarField& x
) const
{
    return amplitude*cos(angle(phase, t, x));
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::AiryCoeffs::vi
(
    const label i,
    const scalar phase,
    const scalar t,
    const vector2DField& xz
) const
{
    const scalarField x(xz.component(0));
    const scalarField z(xz.component(1));

    const scalarField phi(angle(phase, t, x));

    const scalar kzGreat = log(i*great);
    const scalarField kz(min(max(k()*z, - kzGreat), kzGreat));

    if (deep())
    {
        return i*exp(kz)*zip(cos(i*phi), sin(i*phi));
    }
    else
    {
        const scalar kd = k()*depth;
        const scalarField kdz(max(scalar(0), kd + kz));

        return i*zip(cosh(i*kdz)*cos(i*phi), sinh(i*kdz)*sin(i*phi))/sinh(kd);
    }
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::AiryCoeffs::velocity
(
    const scalar phase,
    const scalar t,
    const vector2DField& xz
) const
{
    const scalar ka = k()*amplitude;

    return celerity()*ka*vi(1, phase, t, xz);
}


// ************************************************************************* //
