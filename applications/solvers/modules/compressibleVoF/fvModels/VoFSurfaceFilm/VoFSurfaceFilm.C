/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "VoFSurfaceFilm.H"
#include "compressibleTwoPhaseMixture.H"
#include "fvmSup.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(VoFSurfaceFilm, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            VoFSurfaceFilm,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFSurfaceFilm::VoFSurfaceFilm
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
        thermoSurfaceFilm::typeName,
        mesh,
        mesh.lookupObject<uniformDimensionedVectorField>("g"),
        surfaceFilm::typeName
    ),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::VoFSurfaceFilm::addSupFields() const
{
    return wordList
    (
        {
            surfaceFilm_.rhoPrimary().name(),
            surfaceFilm_.UPrimary().name(),
            surfaceFilm_.TPrimary().name()
        }
    );
}


Foam::scalar Foam::fv::VoFSurfaceFilm::maxDeltaT() const
{
    return surfaceFilm_.maxDeltaT();
}


void Foam::fv::VoFSurfaceFilm::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    surfaceFilm_.evolve();

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::VoFSurfaceFilm::addSup
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

    if (fieldName == surfaceFilm_.rhoPrimary().name())
    {
        eqn += surfaceFilm_.Srho();
    }
    else if (fieldName == surfaceFilm_.TPrimary().name())
    {
        const volScalarField::Internal Cv(surfaceFilm_.primaryThermo().Cv());

        eqn +=
            surfaceFilm_.Sh()()/Cv
          + surfaceFilm_.Srho()
           *(eqn.psi() - surfaceFilm_.primaryThermo().he()/Cv);
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFSurfaceFilm::addSup
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

    eqn += surfaceFilm_.SU();
}


void Foam::fv::VoFSurfaceFilm::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fv::VoFSurfaceFilm::mapMesh(const polyMeshMap& map)
{
    NotImplemented;
}


void Foam::fv::VoFSurfaceFilm::distribute(const polyDistributionMap&)
{
    NotImplemented;
}


bool Foam::fv::VoFSurfaceFilm::movePoints()
{
    return true;
}


// ************************************************************************* //
