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
#include "uniformDimensionedFields.H"
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
    surfaceFilm_
    (
        regionModels::surfaceFilmModels::thermoSingleLayer::typeName,
        mesh,
        mesh.lookupObject<uniformDimensionedVectorField>("g"),
        "surfaceFilm"
    ),
    fieldNames_
    (
        {
            mesh.foundObject<volScalarField>
            (
                IOobject::groupName("rho", surfaceFilm_.phaseName())
            )
          ? IOobject::groupName("rho", surfaceFilm_.phaseName())
          : surfaceFilm_.primaryThermo().rho()().name(),
            surfaceFilm_.UPrimary().name(),
            surfaceFilm_.primaryThermo().he().name()
        }
    ),
    curTimeIndex_(-1)
{
    if (isA<basicSpecieMixture>(surfaceFilm_.primaryThermo()))
    {
        const basicSpecieMixture& composition =
            refCast<const basicSpecieMixture>(surfaceFilm_.primaryThermo());

        const PtrList<volScalarField>& Y = composition.Y();

        forAll(Y, i)
        {
            if (composition.solve(i))
            {
                fieldNames_.append(Y[i].name());
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::surfaceFilm::addSupFields() const
{
    return fieldNames_;
}


void Foam::fv::surfaceFilm::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    surfaceFilm_.evolve();

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

    if (fieldName == fieldNames_[0])
    {
        eqn += surfaceFilm_.Srho();
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

    if (fieldName == fieldNames_[0])
    {
        eqn += surfaceFilm_.Srho();
    }
    else if (fieldName == fieldNames_[2])
    {
        eqn += surfaceFilm_.Sh();
    }
    else if
    (
        isA<basicSpecieMixture>(surfaceFilm_.primaryThermo())
     && refCast<const basicSpecieMixture>(surfaceFilm_.primaryThermo()).contains
        (
            eqn.psi().name()
        )
    )
    {
        eqn += surfaceFilm_.SYi
        (
            refCast<const basicSpecieMixture>(surfaceFilm_.primaryThermo())
           .index(eqn.psi())
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

    if (fieldName == fieldNames_[1])
    {
        eqn += surfaceFilm_.SU();
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


// ************************************************************************* //
