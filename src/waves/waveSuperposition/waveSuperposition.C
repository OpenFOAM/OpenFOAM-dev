/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "waveSuperposition.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::waveSuperposition::transformation
(
    const vectorField& p,
    tensor& axes,
    scalar& u,
    vectorField& xyz
) const
{
    const uniformDimensionedVectorField& g =
        db_.lookupObject<uniformDimensionedVectorField>("g");
    const scalar magG = mag(g.value());
    const vector gHat = g.value()/magG;

    const vector dSurf = direction_ - gHat*(gHat & direction_);
    const scalar magDSurf = mag(dSurf);
    const vector dSurfHat = direction_/magDSurf;

    axes = tensor(dSurfHat, - gHat ^ dSurfHat, - gHat);

    u = speed_*magDSurf;

    xyz = axes & (p - origin_);
}


Foam::tmp<Foam::scalarField> Foam::waveSuperposition::elevation
(
    const scalar t,
    const vector2DField& xy
) const
{
    scalarField result(xy.size(), 0);

    forAll(waveModels_, wavei)
    {
        const vector2D d(cos(waveAngles_[wavei]), sin(waveAngles_[wavei]));
        result += waveModels_[wavei].elevation(t, d.x()*speed_, d & xy);
    }

    return scale(xy)*result;
}


Foam::tmp<Foam::vectorField> Foam::waveSuperposition::velocity
(
    const scalar t,
    const vectorField& xyz
) const
{
    vectorField result(xyz.size(), vector::zero);

    forAll(waveModels_, wavei)
    {
        const vector2D d(cos(waveAngles_[wavei]), sin(waveAngles_[wavei]));
        const vector2DField xz
        (
            zip
            (
                d & zip(xyz.component(0), xyz.component(1)),
                tmp<scalarField>(xyz.component(2))
            )
        );
        const vector2DField uw
        (
            waveModels_[wavei].velocity(t, d.x()*speed_, xz)
        );
        result += zip
        (
            d.x()*uw.component(0),
            d.y()*uw.component(0),
            uw.component(1)
        );
    }

    tmp<scalarField> s = scale(zip(xyz.component(0), xyz.component(1)));

    return s*result;
}


Foam::tmp<Foam::scalarField> Foam::waveSuperposition::pressure
(
    const scalar t,
    const vectorField& xyz
) const
{
    scalarField result(xyz.size(), 0);

    forAll(waveModels_, wavei)
    {
        const vector2D d(cos(waveAngles_[wavei]), sin(waveAngles_[wavei]));
        const vector2DField xz
        (
            zip
            (
                d & zip(xyz.component(0), xyz.component(1)),
                tmp<scalarField>(xyz.component(2))
            )
        );
        const vector2DField uw
        (
            waveModels_[wavei].velocity(t, d.x()*speed_, xz)
        );
        result += waveModels_[wavei].pressure(t, d.x()*speed_, xz);
    }

    tmp<scalarField> s = scale(zip(xyz.component(0), xyz.component(1)));

    return s*result;
}


Foam::tmp<Foam::scalarField> Foam::waveSuperposition::scale
(
    const vector2DField& xy
) const
{
    tmp<scalarField> tResult(new scalarField(xy.size(), 1));
    scalarField& result = tResult.ref();

    if (scale_.valid())
    {
        const scalarField x(xy.component(0));
        forAll(result, i)
        {
            result[i] *= scale_->value(x[i]);
        }
    }

    if (crossScale_.valid())
    {
        const scalarField y(xy.component(1));
        forAll(result, i)
        {
            result[i] *= crossScale_->value(y[i]);
        }
    }

    return tResult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSuperposition::waveSuperposition(const objectRegistry& db)
:
    db_(db),
    origin_(vector::zero),
    direction_(vector(1, 0, 0)),
    speed_(0),
    waveModels_(),
    waveAngles_(),
    ramp_(),
    scale_(),
    crossScale_(),
    heightAboveWave_(false)
{}


Foam::waveSuperposition::waveSuperposition(const waveSuperposition& waves)
:
    db_(waves.db_),
    origin_(waves.origin_),
    direction_(waves.direction_),
    speed_(waves.speed_),
    waveModels_(waves.waveModels_),
    waveAngles_(waves.waveAngles_),
    ramp_(waves.ramp_, false),
    scale_(waves.scale_, false),
    crossScale_(waves.crossScale_, false),
    heightAboveWave_(waves.heightAboveWave_)
{}


Foam::waveSuperposition::waveSuperposition
(
    const objectRegistry& db,
    const dictionary& dict
)
:
    db_(db),
    origin_(dict.lookup("origin")),
    direction_(dict.lookup("direction")),
    speed_(readScalar(dict.lookup("speed"))),
    waveModels_(),
    waveAngles_(),
    ramp_
    (
        dict.found("ramp")
      ? Function1<scalar>::New("ramp", dict)
      : autoPtr<Function1<scalar>>()
    ),
    scale_
    (
        dict.found("scale")
      ? Function1<scalar>::New("scale", dict)
      : autoPtr<Function1<scalar>>()
    ),
    crossScale_
    (
        dict.found("crossScale")
      ? Function1<scalar>::New("crossScale", dict)
      : autoPtr<Function1<scalar>>()
    ),
    heightAboveWave_(dict.lookupOrDefault<Switch>("heightAboveWave", false))
{
    const PtrList<entry> waveEntries(dict.lookup("waves"));

    waveModels_.setSize(waveEntries.size());
    waveAngles_.setSize(waveEntries.size());

    forAll(waveEntries, wavei)
    {
        const dictionary waveDict = waveEntries[wavei].dict();
        waveModels_.set
        (
            wavei,
            waveModel::New(waveDict.dictName(), db, waveDict)
        );
        waveAngles_[wavei] = readScalar(waveDict.lookup("angle"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveSuperposition::~waveSuperposition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveSuperposition::height
(
    const scalar t,
    const vectorField& p
) const
{
    tensor axes;
    scalar u;
    vectorField xyz(p.size());
    transformation(p, axes, u, xyz);

    return
        xyz.component(2)
      - elevation(t, zip(xyz.component(0), xyz.component(1)));
}


Foam::tmp<Foam::vectorField> Foam::waveSuperposition::ULiquid
(
    const scalar t,
    const vectorField& p
) const
{
    tensor axes;
    scalar u;
    vectorField xyz(p.size());
    transformation(p, axes, u, xyz);

    if (heightAboveWave_)
    {
        xyz.replace(2, height(t, p));
    }

    return UMean(t) + (velocity(t, xyz) & axes);
}


Foam::tmp<Foam::vectorField> Foam::waveSuperposition::UGas
(
    const scalar t,
    const vectorField& p
) const
{
    tensor axes;
    scalar u;
    vectorField xyz(p.size());
    transformation(p, axes, u, xyz);

    axes = tensor(- axes.x(), - axes.y(), axes.z());

    if (heightAboveWave_)
    {
        xyz.replace(2, height(t, p));
    }

    return UMean(t) + (velocity(t, xyz) & axes);
}


Foam::tmp<Foam::scalarField> Foam::waveSuperposition::pLiquid
(
    const scalar t,
    const vectorField& p
) const
{
    tensor axes;
    scalar u;
    vectorField xyz(p.size());
    transformation(p, axes, u, xyz);

    if (heightAboveWave_)
    {
        xyz.replace(2, height(t, p));
    }

    return pressure(t, xyz);
}


Foam::tmp<Foam::scalarField> Foam::waveSuperposition::pGas
(
    const scalar t,
    const vectorField& p
) const
{
    return - pLiquid(t, p);
}


void Foam::waveSuperposition::write(Ostream& os) const
{
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction") << direction_ << token::END_STATEMENT << nl;
    os.writeKeyword("speed") << speed_ << token::END_STATEMENT << nl;
    os.writeKeyword("waves") << nl << token::BEGIN_LIST << nl << incrIndent;
    forAll(waveModels_, wavei)
    {
        os.writeKeyword(waveModels_[wavei].type()) << nl << indent
            << token::BEGIN_BLOCK << nl << incrIndent;
        waveModels_[wavei].write(os);
        os.writeKeyword("angle") << waveAngles_[wavei] << token::END_STATEMENT
            << nl << decrIndent << indent << token::END_BLOCK << nl;
    }
    os  << decrIndent << token::END_LIST << token::END_STATEMENT << nl;
    if (ramp_.valid())
    {
        ramp_->writeData(os);
    }
    if (scale_.valid())
    {
        scale_->writeData(os);
    }
    if (crossScale_.valid())
    {
        crossScale_->writeData(os);
    }
    if (heightAboveWave_)
    {
        os.writeKeyword("heightAboveWave") << heightAboveWave_
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
