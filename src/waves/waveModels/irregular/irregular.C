/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "irregular.H"
#include "mathematicalConstants.H"
#include "OSspecific.H"
#include "randomGenerator.H"
#include "rawSetWriter.H"
#include "writeFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(irregular, 0);
    addToRunTimeSelectionTable(waveModel, irregular, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::waveModels::AiryCoeffs
Foam::waveModels::irregular::coeffs(const label i) const
{
    return AiryCoeffs(depth_, amplitudes_[i], lengths_[i], g());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::irregular::irregular(const irregular& wave)
:
    waveModel(wave),
    depth_(wave.depth_),
    spectrum_(wave.spectrum_, false),
    n_(wave.n_),
    span_(wave.span_),
    formatter_(wave.formatter_, false),
    amplitudes_(wave.amplitudes_),
    lengths_(wave.lengths_),
    phases_(wave.phases_)
{}


Foam::waveModels::irregular::irregular
(
    const dictionary& dict,
    const scalar g,
    const word& modelName
)
:
    waveModel(dict, g),
    depth_(dict.lookupOrDefault<scalar>("depth", great)),
    spectrum_(waveSpectrum::New(dict, g)),
    n_(dict.lookup<label>("n")),
    span_(dict.lookupOrDefault("span", Pair<scalar>(0.01, 0.99))),
    formatter_
    (
        dict.found("setFormat")
      ? setWriter::New(dict.lookup("setFormat"), dict).ptr()
      : nullptr
    ),
    amplitudes_(n_),
    lengths_(n_),
    phases_(n_)
{
    // Get the bounds of the frequency space
    const Pair<scalar> f01
    (
        spectrum_->fFraction(span_.first()),
        spectrum_->fFraction(span_.second())
    );

    // Create a discretisation in frequency space
    const scalarField x(scalarField(scalarList(identityMap(n_ + 1)))/n_);
    const scalarField f((1 - x)*f01.first() + x*f01.second());

    // Evaluate the integrals of the distribution
    const scalarField integralS(spectrum_->integralS(f));
    const scalarField integralFS(spectrum_->integralFS(f));

    // Set the wave parameters using means constructed from the integrals over
    // the frequency intervals
    scalarList periods(n_);
    forAll(amplitudes_, i)
    {
        const label j0 = n_ - 1 - i, j1 = n_ - i;

        // Integrating the of power spectral density gives the square of the
        // root-mean-squared value. So, a factor of two is needed to convert to
        // peak amplitude.
        amplitudes_[i] = sqrt(2*(integralS[j1] - integralS[j0]));

        // Frequency is the integrated average, Integral(f*S)/Integral(S), for
        // this interval. Period is the reciprocal of this value.
        periods[i] =
            (integralS[j1] - integralS[j0])/(integralFS[j1] - integralFS[j0]);

        // Let the Airy calculation engine convert the above to wavelength
        lengths_[i] =
            AiryCoeffs
            (
                depth_,
                amplitudes_[i],
                periods[i],
                g,
                AiryCoeffs::celerity
            ).length;
    }

    // Set random phases
    randomGenerator rndGen(0, true);
    forAll(phases_, i)
    {
        phases_[i] = rndGen.scalarAB(0, constant::mathematical::twoPi);
    }

    // Report
    const scalar avgAmplitude =
        sqrt(2*(integralS.last() - integralS.first())/(f.last() - f.first()));
    const scalar avgPeroid =
        (integralS.last() - integralS.first())
       /(integralFS.last() - integralFS.first());
    const scalar avgLength =
        AiryCoeffs(depth_, avgAmplitude, avgPeroid, g, AiryCoeffs::celerity)
       .length;
    const string::size_type pad =
        waveModel::typeName.length() + modelName.length() + 2;
    Info<< waveModel::typeName << ": " << modelName
        << ": min/average/max period = "
        << periods.first() << '/' << avgPeroid << '/' << periods.last()
        << endl << string(pad, ' ').c_str()
        << "  min/average/max length = "
        << lengths_.first() << '/' << avgLength << '/' << lengths_.last()
        << endl;

    // Write out the continuous distribution
    if (Pstream::master() && formatter_.valid())
    {
        static const scalar t = 0.1;

        const Pair<scalar> g01
        (
            max(rootVSmall, (1 + t)*(f01.first() - t*f01.second())),
            - t*f01.first() + (1 + t)*f01.second()
        );

        const label nn = 1024;
        const scalarField xx(scalarField(scalarList(identityMap(nn + 1)))/nn);
        const scalarField ff((1 - xx)*g01.first() + xx*g01.second());

        formatter_->write
        (
            cwd()
           /functionObjects::writeFile::outputPrefix
           /waveModel::typeName + "s"
           /modelName + spectrum_->type().capitalise(),
            "continuous",
            coordSet(true, "f", ff),
            "S",
            spectrum_->S(ff)()
        );
    }

    // Write out the sampled distribution
    if (Pstream::master() && formatter_.valid())
    {
        scalarField ff(3*n_ + 1), SS(3*n_ + 1);

        forAll(f, i)
        {
            if (i != 0)
            {
                ff[3*i - 1] = f[i];
                SS[3*i - 1] = sqr(amplitudes_[n_ - i])/2/(f[i] - f[i - 1]);
            }

            ff[3*i] = f[i];
            SS[3*i] = 0;

            if (i != n_)
            {
                ff[3*i + 1] = f[i];
                SS[3*i + 1] = sqr(amplitudes_[n_ - 1 - i])/2/(f[i + 1] - f[i]);
            }
        }

        formatter_->write
        (
            cwd()
           /functionObjects::writeFile::outputPrefix
           /waveModel::typeName + "s"
           /modelName + spectrum_->type().capitalise(),
            "sampled",
            coordSet(true, "f", ff),
            "S",
            SS
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::irregular::~irregular()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::irregular::celerity() const
{
    return coeffs(amplitudes_.size() - 1).celerity();
}


Foam::tmp<Foam::scalarField> Foam::waveModels::irregular::elevation
(
    const scalar t,
    const scalarField& x
) const
{
    tmp<scalarField> tResult(new scalarField(x.size(), scalar(0)));
    scalarField& result = tResult.ref();

    forAll(amplitudes_, i)
    {
        result += coeffs(i).elevation(phases_[i], t, x);
    }

    return tResult;
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::irregular::velocity
(
    const scalar t,
    const vector2DField& xz
) const
{
    tmp<vector2DField> tResult(new vector2DField(xz.size(), vector2D::zero));
    vector2DField& result = tResult.ref();

    forAll(amplitudes_, i)
    {
        result += coeffs(i).velocity(phases_[i], t, xz);
    }

    return tResult;
}


void Foam::waveModels::irregular::write(Ostream& os) const
{
    waveModel::write(os);

    if (!coeffs(amplitudes_.size() - 1).deep())
    {
        writeEntry(os, "depth", depth_);
    }

    writeEntry(os, "spectrum", spectrum_->type());

    os  << indent << word(spectrum_->type() + "Coeffs") << nl
        << indent << token::BEGIN_BLOCK << nl << incrIndent;
    spectrum_->write(os);
    os  << decrIndent
        << indent << token::END_BLOCK << nl;

    writeEntry(os, "n", n_);
    writeEntryIfDifferent(os, "span", Pair<scalar>(0.01, 0.99), span_);
}


// ************************************************************************* //
