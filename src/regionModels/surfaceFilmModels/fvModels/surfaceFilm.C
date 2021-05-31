/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "surfaceFilm.H"
#include "basicSpecieMixture.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(surfaceFilm, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            surfaceFilm,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::surfaceFilm::surfaceFilm
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(sourceName, modelType, dict, mesh),
    primaryThermo_
    (
        mesh.lookupObject<fluidThermo>(basicThermo::dictName)
    ),
    surfaceFilm_
    (
        regionModels::surfaceFilmModel::New
        (
            mesh,
            mesh.lookupObject<uniformDimensionedVectorField>("g")
        )
    ),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::surfaceFilm::addSupFields() const
{
    wordList fieldNames({"rho", "U", primaryThermo_.he().name()});

    if (isA<basicSpecieMixture>(primaryThermo_))
    {
        const basicSpecieMixture& composition =
            refCast<const basicSpecieMixture>(primaryThermo_);

        const PtrList<volScalarField>& Y = composition.Y();

        forAll(Y, i)
        {
            if (composition.solve(i))
            {
                fieldNames.append(Y[i].name());
            }
        }
    }

    return fieldNames;
}


void Foam::fv::surfaceFilm::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    surfaceFilm_->evolve();

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::surfaceFilm::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == "rho")
    {
        eqn += surfaceFilm_->Srho();
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::surfaceFilm::addSup
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

    if (fieldName == "rho")
    {
        eqn += surfaceFilm_->Srho();
    }
    else if (fieldName == primaryThermo_.he().name())
    {
        eqn += surfaceFilm_->Sh();
    }
    else if
    (
        isA<basicSpecieMixture>(primaryThermo_)
     && refCast<const basicSpecieMixture>(primaryThermo_).contains
        (
            eqn.psi().name()
        )
    )
    {
        eqn += surfaceFilm_->SYi
        (
            refCast<const basicSpecieMixture>(primaryThermo_).index(eqn.psi())
        );
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::surfaceFilm::addSup
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

    if (fieldName == "U")
    {
        eqn += surfaceFilm_->SU();
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


// ************************************************************************* //
