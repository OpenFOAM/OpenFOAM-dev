/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "Airy.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Airy, 0);
    addToRunTimeSelectionTable(waveModel, Airy, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Airy::length
(
    const dictionary& dict,
    const scalar depth,
    const scalar amplitude,
    const scalar g,
    scalar (*modelCelerity)(scalar, scalar, scalar, scalar)
)
{
    const bool haveLength = dict.found("length");
    const bool havePeriod = dict.found("period");

    if (haveLength == havePeriod)
    {
        FatalIOErrorInFunction(dict)
            << "Exactly one of either length or period must be specified"
            << exit(FatalIOError);
    }

    if (haveLength)
    {
        return dict.lookup<scalar>("length");
    }
    else
    {
        const scalar period = dict.lookup<scalar>("period");

        scalar length0 = 0;
        scalar length1 = 1.5*g*sqr(period)/(2*constant::mathematical::pi);

        // Bisect down to round off error to find the wave length
        for (label i = 0; i < ceil(std::log2(1/small)); ++ i)
        {
            const scalar length = (length0 + length1)/2;

            const scalar t = length/modelCelerity(depth, amplitude, length, g);

            (t < period ? length0 : length1) = length;
        }

        return (length0 + length1)/2;
    }
}


// * * * * * * * * * * Static Protected Member Functions  * * * * * * ** * * //

Foam::scalar Foam::waveModels::Airy::k(const scalar length)
{
    return 2*Foam::constant::mathematical::pi/length;
}


bool Foam::waveModels::Airy::deep(const scalar depth, const scalar length)
{
    return depth*k(length) > log(great);
}


Foam::scalar Foam::waveModels::Airy::celerity
(
    const scalar depth,
    const scalar amplitude,
    const scalar length,
    const scalar g
)
{
    return sqrt(g/k(length)*tanh(k(length)*depth));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::Airy::angle
(
    const scalar t,
    const scalarField& x
) const
{
    return phase_ + k()*(x - celerity()*t);
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Airy::vi
(
    const label i,
    const scalar t,
    const vector2DField& xz
) const
{
    const scalarField x(xz.component(0));
    const scalarField z(xz.component(1));

    const scalarField phi(angle(t, x));

    const scalar kzGreat = log(i*great);
    const scalarField kz(min(max(k()*z, - kzGreat), kzGreat));

    if (deep())
    {
        return i*exp(kz)*zip(cos(i*phi), sin(i*phi));
    }
    else
    {
        const scalar kd = k()*depth();
        const scalarField kdz(max(scalar(0), kd + kz));

        return i*zip(cosh(i*kdz)*cos(i*phi), sinh(i*kdz)*sin(i*phi))/sinh(kd);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Airy::Airy(const Airy& wave)
:
    waveModel(wave),
    depth_(wave.depth_),
    amplitude_(wave.amplitude_, false),
    length_(wave.length_),
    phase_(wave.phase_)
{}


Foam::waveModels::Airy::Airy
(
    const dictionary& dict,
    const scalar g,
    const word& modelName,
    scalar (*modelCelerity)(scalar, scalar, scalar, scalar)
)
:
    waveModel(dict, g),
    depth_(dict.lookupOrDefault<scalar>("depth", great)),
    amplitude_(Function1<scalar>::New("amplitude", dict)),
    length_(length(dict, depth_, amplitude(), g, modelCelerity)),
    phase_(dict.lookup<scalar>("phase"))
{
    const scalar c = modelCelerity(depth(), amplitude(), length(), g);

    Info<< waveModel::typeName << ": " << modelName
        << ": period = " << length_/c
        << ", length = " << length_ << endl;;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::Airy::~Airy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::Airy::elevation
(
    const scalar t,
    const scalarField& x
) const
{
    return amplitude(t)*cos(angle(t, x));
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Airy::velocity
(
    const scalar t,
    const vector2DField& xz
) const
{
    const scalar ka = k()*amplitude(t);

    return Airy::celerity()*ka*vi(1, t, xz);
}


void Foam::waveModels::Airy::write(Ostream& os) const
{
    waveModel::write(os);

    if (!deep())
    {
        writeEntry(os, "depth", depth_);
    }
    writeEntry(os, amplitude_());
    writeEntry(os, "length", length_);
    writeEntry(os, "phase", phase_);
}


// ************************************************************************* //
