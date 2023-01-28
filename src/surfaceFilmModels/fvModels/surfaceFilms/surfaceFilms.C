/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "surfaceFilms.H"
#include "uniformDimensionedFields.H"
#include "basicSpecieMixture.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(surfaceFilms, 0);

        addToRunTimeSelectionTable(fvModel, surfaceFilms, dictionary);

        // !!! Backwards compatible lookup name
        addNamedToRunTimeSelectionTable
        (
            fvModel,
            surfaceFilms,
            dictionary,
            surfaceFilm
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::surfaceFilms::surfaceFilms
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    surfaceFilms_(),
    surfaceFilmPrimaryRhoNames_(),
    fieldNames_(),
    curTimeIndex_(-1)
{
    const wordList surfaceFilmNames =
        dict.lookupOrDefault<wordList>
        (
            "surfaceFilms",
            wordList(1, surfaceFilm::typeName)
        );

    surfaceFilms_.resize(surfaceFilmNames.size());
    surfaceFilmPrimaryRhoNames_.resize(surfaceFilmNames.size());

    wordHashSet fieldNamesSet;

    forAll(surfaceFilmNames, i)
    {
        surfaceFilms_.set
        (
            i,
            new thermoSurfaceFilm
            (
                thermoSurfaceFilm::typeName,
                mesh,
                mesh.lookupObject<uniformDimensionedVectorField>("g"),
                surfaceFilmNames[i]
            )
        );

        const basicThermo& primaryThermo = surfaceFilms_[i].primaryThermo();

        const word transportRhoName =
            IOobject::groupName("rho", surfaceFilms_[i].phaseName());

        surfaceFilmPrimaryRhoNames_[i] =
            mesh.foundObject<volScalarField>(transportRhoName)
          ? transportRhoName
          : primaryThermo.rho()().name();

        fieldNamesSet.insert(surfaceFilmPrimaryRhoNames_[i]);
        fieldNamesSet.insert(surfaceFilms_[i].UPrimary().name());
        fieldNamesSet.insert(primaryThermo.he().name());

        if (isA<basicSpecieMixture>(primaryThermo))
        {
            const basicSpecieMixture& mixture =
                refCast<const basicSpecieMixture>(primaryThermo);

            forAll(mixture.Y(), Yi)
            {
                if (mixture.solve(Yi))
                {
                    fieldNamesSet.insert(mixture.Y()[Yi].name());
                }
            }
        }
    }

    fieldNames_ = fieldNamesSet.toc();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::surfaceFilms::addSupFields() const
{
    return fieldNames_;
}


Foam::scalar Foam::fv::surfaceFilms::maxDeltaT() const
{
    scalar result = vGreat;
    forAll(surfaceFilms_, i)
    {
        result = min(result, surfaceFilms_[i].maxDeltaT());
    }
    return result;
}


void Foam::fv::surfaceFilms::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    forAll(surfaceFilms_, i)
    {
        surfaceFilms_[i].evolve();
    }

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::surfaceFilms::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    forAll(surfaceFilms_, i)
    {
        if (eqn.psi().group() == surfaceFilms_[i].phaseName())
        {
            if (fieldName == surfaceFilmPrimaryRhoNames_[i])
            {
                eqn += surfaceFilms_[i].Srho();
            }
            else
            {
                FatalErrorInFunction
                    << "Support for field " << fieldName
                    << " is not implemented" << exit(FatalError);
            }
        }
    }
}


void Foam::fv::surfaceFilms::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    forAll(surfaceFilms_, i)
    {
        const basicThermo& primaryThermo = surfaceFilms_[i].primaryThermo();

        if (eqn.psi().group() == surfaceFilms_[i].phaseName())
        {
            if (fieldName == surfaceFilmPrimaryRhoNames_[i])
            {
                forAll(surfaceFilms_, i)
                {
                    eqn += surfaceFilms_[i].Srho();
                }
            }
            else if (fieldName == primaryThermo.he().name())
            {
                forAll(surfaceFilms_, i)
                {
                    eqn += surfaceFilms_[i].Sh();
                }
            }
            else if
            (
                isA<basicSpecieMixture>(primaryThermo)
             && refCast<const basicSpecieMixture>(primaryThermo).contains
                (
                    eqn.psi().name()
                )
            )
            {
                eqn += surfaceFilms_[i].SYi
                (
                    refCast<const basicSpecieMixture>(primaryThermo).index
                    (
                        eqn.psi()
                    )
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Support for field " << fieldName
                    << " is not implemented" << exit(FatalError);
            }
        }
    }
}


void Foam::fv::surfaceFilms::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    forAll(surfaceFilms_, i)
    {
        if (eqn.psi().group() == surfaceFilms_[i].phaseName())
        {
            if (fieldName == surfaceFilms_[i].UPrimary().name())
            {
                eqn += surfaceFilms_[i].SU();
            }
            else
            {
                FatalErrorInFunction
                    << "Support for field " << fieldName
                    << " is not implemented" << exit(FatalError);
            }
        }
    }
}


void Foam::fv::surfaceFilms::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fv::surfaceFilms::mapMesh(const polyMeshMap& map)
{
    NotImplemented;
}


void Foam::fv::surfaceFilms::distribute(const polyDistributionMap&)
{
    NotImplemented;
}


bool Foam::fv::surfaceFilms::movePoints()
{
    return true;
}


// ************************************************************************* //
