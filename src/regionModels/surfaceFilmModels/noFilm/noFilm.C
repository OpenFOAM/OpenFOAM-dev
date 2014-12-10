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

#include "noFilm.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noFilm, 0);
addToRunTimeSelectionTable(surfaceFilmModel, noFilm, mesh);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool noFilm::read()
{
    if (surfaceFilmModel::read())
    {
        // no additional info to read
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noFilm::noFilm
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType
)
:
    surfaceFilmModel(modelType, mesh, g, regionType)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noFilm::~noFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void noFilm::addSources
(
    const label,
    const label,
    const scalar,
    const vector&,
    const scalar,
    const scalar
)
{
    // do nothing
}


const volScalarField& noFilm::delta() const
{
    FatalErrorIn("const volScalarField& noFilm::delta() const")
        << "delta field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::alpha() const
{
    FatalErrorIn("const volScalarField& noFilm::alpha() const")
        << "alpha field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volVectorField& noFilm::U() const
{
    FatalErrorIn("const volVectorField& noFilm::U() const")
        << "U field not available for " << type() << abort(FatalError);

    return volVectorField::null();
}


const volVectorField& noFilm::Us() const
{
    FatalErrorIn("const volVectorField& noFilm::Us() const")
        << "Us field not available for " << type() << abort(FatalError);

    return volVectorField::null();
}


const volVectorField& noFilm::Uw() const
{
    FatalErrorIn("const volVectorField& noFilm::Uw() const")
        << "Uw field not available for " << type() << abort(FatalError);

    return volVectorField::null();
}


const volScalarField& noFilm::rho() const
{
    FatalErrorIn("const volScalarField& noFilm::rho() const")
        << "rho field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::T() const
{
    FatalErrorIn("const volScalarField& noFilm::T() const")
        << "T field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::Ts() const
{
    FatalErrorIn("const volScalarField& noFilm::Ts() const")
        << "Ts field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::Tw() const
{
    FatalErrorIn("const volScalarField& noFilm::Tw() const")
        << "Tw field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::Cp() const
{
    FatalErrorIn("const volScalarField& noFilm::Cp() const")
        << "Cp field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::kappa() const
{
    FatalErrorIn("const volScalarField& noFilm::kappa() const")
        << "kappa field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::sigma() const
{
    FatalErrorIn("const volScalarField& noFilm::sigma() const")
        << "sigma field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


tmp<volScalarField> noFilm::primaryMassTrans() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "noFilm::primaryMassTrans",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );
}


const volScalarField& noFilm::cloudMassTrans() const
{
    FatalErrorIn("const volScalarField& noFilm::cloudMassTrans() const")
        << "cloudMassTrans field not available for " << type()
        << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& noFilm::cloudDiameterTrans() const
{
    FatalErrorIn("const volScalarField& noFilm::cloudDiameterTrans() const")
        << "cloudDiameterTrans field not available for " << type()
        << abort(FatalError);

    return volScalarField::null();
}


tmp<DimensionedField<scalar, volMesh> > noFilm::Srho() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "noFilm::Srho",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );
}


tmp<DimensionedField<scalar, volMesh> > noFilm::Srho(const label i) const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "noFilm::Srho(" + Foam::name(i) + ")",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );
}


tmp<DimensionedField<scalar, volMesh> > noFilm::Sh() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "noFilm::Sh",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
