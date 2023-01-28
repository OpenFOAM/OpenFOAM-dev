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

#include "forcing.H"
#include "fvMatrix.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(forcing, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::forcing::readCoeffs()
{
    lambda_ =
        dimensionedScalar
        (
            lambda_.name(),
            lambda_.dimensions(),
            coeffs().lookup(lambda_.name())
        );

    const bool foundScale = coeffs().found("scale");
    const bool foundOgn = coeffs().found("origin");
    const bool foundDir = coeffs().found("direction");
    const bool foundOgns = coeffs().found("origins");
    const bool foundDirs = coeffs().found("directions");
    const bool foundAll =
        foundScale
     && (
            (foundOgn && foundDir && !foundOgns && !foundDirs)
         || (!foundOgn && !foundDir && foundOgns && foundDirs)
        );
     const bool foundAny =
        foundScale || foundOgn || foundDir || foundOgns || foundDirs;

    if (!foundAll)
    {
        scale_ = autoPtr<Function1<scalar>>();
        origins_.clear();
        directions_.clear();
    }

    if (foundAll)
    {
        scale_ = Function1<scalar>::New("scale", coeffs());
        if (foundOgn)
        {
            origins_.setSize(1);
            directions_.setSize(1);
            coeffs().lookup("origin") >> origins_.last();
            coeffs().lookup("direction") >> directions_.last();
        }
        else
        {
            coeffs().lookup("origins") >> origins_;
            coeffs().lookup("directions") >> directions_;

            if
            (
                origins_.size() == 0
             || directions_.size() == 0
             || origins_.size() != directions_.size()
            )
            {
                FatalErrorInFunction
                    << "The same, non-zero number of origins and "
                    << "directions must be provided" << exit(FatalError);
            }
        }
        forAll(directions_, i)
        {
            directions_[i] /= mag(directions_[i]);
        }
    }

    if (!foundAll && foundAny)
    {
        WarningInFunction
            << "The scaling specification is incomplete. \"scale\", "
            << "\"origin\" and \"direction\" (or \"origins\" and "
            << "\"directions\"), must all be specified in order to scale "
            << "the forcing. The forcing will be applied uniformly across "
            << "the cell set." << endl << endl;
    }
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::forcing::scale() const
{
    tmp<volScalarField::Internal> tscale
    (
        volScalarField::Internal::New
        (
            typedName("scale"),
            mesh(),
            dimensionedScalar(dimless, scale_.valid() ? 0 : 1)
        )
    );

    scalarField& scale = tscale.ref();

    forAll(origins_, i)
    {
        const vectorField& c = mesh().cellCentres();
        const scalarField x((c - origins_[i]) & directions_[i]);
        scale = max(scale, scale_->value(x));
    }

    // Write out the force coefficient for debugging
    if (debug && mesh().time().writeTime())
    {
        tscale->write();
    }

    return tscale;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::forcing::forceCoeff
(
    const volScalarField::Internal& scale
) const
{
    tmp<volScalarField::Internal> tforceCoeff
    (
        volScalarField::Internal::New(typedName("forceCoeff"), lambda_*scale)
    );

    // Write out the force coefficient for debugging
    if (debug && mesh().time().writeTime())
    {
        tforceCoeff->write();
    }

    return tforceCoeff;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::forcing::forceCoeff() const
{
    return forceCoeff(scale());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::forcing::forcing
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    lambda_("lambda", dimless/dimTime, NaN),
    scale_(nullptr),
    origins_(),
    directions_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::forcing::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
