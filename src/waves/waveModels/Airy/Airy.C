/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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

Foam::scalar Foam::waveModels::Airy::readLength
(
    const dictionary& dict,
    const scalar depth,
    const scalar amplitude,
    const scalar g,
    scalar (*celerityPtr)(const AiryCoeffs&)
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
        return
            AiryCoeffs
            (
                depth,
                amplitude,
                dict.lookup<scalar>("period"),
                g,
                celerityPtr
            ).length;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Airy::celerity(const AiryCoeffs& coeffs)
{
    return coeffs.celerity();
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
    scalar (*celerityPtr)(const AiryCoeffs&)
)
:
    waveModel(dict, g),
    depth_(dict.lookupOrDefault<scalar>("depth", great)),
    amplitude_(Function1<scalar>::New("amplitude", dimTime, dimLength, dict)),
    length_(readLength(dict, depth_, amplitude(), g, celerityPtr)),
    phase_(dict.lookup<scalar>("phase"))
{
    const scalar c = celerityPtr(coeffs());

    Info<< waveModel::typeName << ": " << modelName
        << ": period = " << length_/c
        << ", length = " << length_ << endl;;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::Airy::~Airy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Airy::celerity() const
{
    return celerity(coeffs());
}


Foam::tmp<Foam::scalarField> Foam::waveModels::Airy::elevation
(
    const scalar t,
    const scalarField& x
) const
{
    return amplitude(t)*cos(coeffs().angle(phase_, t, x));
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Airy::velocity
(
    const scalar t,
    const vector2DField& xz
) const
{
    const scalar ka = coeffs().k()*amplitude(t);

    return Airy::celerity()*ka*coeffs().vi(1, phase_, t, xz);
}


void Foam::waveModels::Airy::write(Ostream& os) const
{
    waveModel::write(os);

    if (!coeffs().deep())
    {
        writeEntry(os, "depth", depth_);
    }
    writeEntry(os, amplitude_());
    writeEntry(os, "length", length_);
    writeEntry(os, "phase", phase_);
}


// ************************************************************************* //
