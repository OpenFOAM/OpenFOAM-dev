/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
#include "Time.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::waveSuperposition::dictName("waveProperties");

namespace Foam
{
    defineTypeNameAndDebug(waveSuperposition, 0);
    defineRunTimeSelectionTable(waveSuperposition, objectRegistry);
    addToRunTimeSelectionTable
    (
        waveSuperposition,
        waveSuperposition,
        objectRegistry
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::waveSuperposition::transformation
(
    const scalar t,
    const vectorField& p,
    tensor& axes,
    vectorField& xyz
) const
{
    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");
    const scalar magG = mag(g.value());
    const vector gHat = g.value()/magG;

    const vector dSurf = direction_ - gHat*(gHat & direction_);
    const scalar magDSurf = mag(dSurf);
    const vector dSurfHat = direction_/magDSurf;

    axes = tensor(dSurfHat, - gHat ^ dSurfHat, - gHat);

    xyz = axes & (p - origin_ - UMean_->integral(0, t));
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
        result += waveModels_[wavei].elevation(t, d & xy);
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
            waveModels_[wavei].velocity(t, xz)
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
            waveModels_[wavei].velocity(t, xz)
        );
        result += waveModels_[wavei].pressure(t, xz);
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
    IOdictionary
    (
        IOobject
        (
            dictName,
            db.time().constant(),
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    origin_(lookup("origin")),
    direction_(lookup("direction")),
    waveModels_(),
    waveAngles_(),
    UMean_(Function1<vector>::New("UMean", *this)),
    scale_
    (
        found("scale")
      ? Function1<scalar>::New("scale", *this)
      : autoPtr<Function1<scalar>>()
    ),
    crossScale_
    (
        found("crossScale")
      ? Function1<scalar>::New("crossScale", *this)
      : autoPtr<Function1<scalar>>()
    ),
    heightAboveWave_(lookupOrDefault<Switch>("heightAboveWave", false))
{
    const PtrList<entry> waveEntries(lookup("waves"));

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
        waveAngles_[wavei] = waveDict.lookup<scalar>("angle");
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
    vectorField xyz(p.size());
    transformation(t, p, axes, xyz);

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
    vectorField xyz(p.size());
    transformation(t, p, axes, xyz);

    if (heightAboveWave_)
    {
        xyz.replace(2, height(t, p));
    }

    return UMean_->value(t) + (velocity(t, xyz) & axes);
}


Foam::tmp<Foam::vectorField> Foam::waveSuperposition::UGas
(
    const scalar t,
    const vectorField& p
) const
{
    tensor axes;
    vectorField xyz(p.size());
    transformation(t, p, axes, xyz);

    axes = tensor(- axes.x(), - axes.y(), axes.z());

    if (heightAboveWave_)
    {
        xyz.replace(2, height(t, p));
    }

    xyz.replace(2, - xyz.component(2));

    return UMean_->value(t) + (velocity(t, xyz) & axes);
}


Foam::tmp<Foam::scalarField> Foam::waveSuperposition::pLiquid
(
    const scalar t,
    const vectorField& p
) const
{
    tensor axes;
    vectorField xyz(p.size());
    transformation(t, p, axes, xyz);

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
    tensor axes;
    vectorField xyz(p.size());
    transformation(t, p, axes, xyz);

    axes = tensor(- axes.x(), - axes.y(), axes.z());

    if (heightAboveWave_)
    {
        xyz.replace(2, height(t, p));
    }

    xyz.replace(2, - xyz.component(2));

    return pressure(t, xyz);
}


void Foam::waveSuperposition::write(Ostream& os) const
{
    writeEntry(os, "origin", origin_);
    writeEntry(os, "direction", direction_);
    writeKeyword(os, "waves") << nl << token::BEGIN_LIST << nl << incrIndent;
    forAll(waveModels_, wavei)
    {
        writeKeyword(os, waveModels_[wavei].type()) << nl << indent
            << token::BEGIN_BLOCK << nl << incrIndent;
        waveModels_[wavei].write(os);
        writeKeyword(os, "angle") << waveAngles_[wavei] << token::END_STATEMENT
            << nl << decrIndent << indent << token::END_BLOCK << nl;
    }
    os  << decrIndent << token::END_LIST << token::END_STATEMENT << nl;
    writeEntry(os, UMean_());
    if (scale_.valid())
    {
        writeEntry(os, scale_());
    }
    if (crossScale_.valid())
    {
        writeEntry(os, crossScale_());
    }
    if (heightAboveWave_)
    {
        writeEntry(os, "heightAboveWave", heightAboveWave_);
    }
}


// ************************************************************************* //
