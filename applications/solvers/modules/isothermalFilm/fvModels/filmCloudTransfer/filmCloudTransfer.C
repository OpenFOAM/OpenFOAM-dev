/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "filmCloudTransfer.H"
#include "mappedPatchBase.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(filmCloudTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            filmCloudTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::filmCloudTransfer::filmCloudTransfer
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    film_(mesh.lookupObject<solvers::isothermalFilm>(solver::typeName)),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::filmCloudTransfer::addSupFields() const
{
    return wordList
    {
        "pi",
        film_.alpha.name(),
        film_.thermo.he().name(),
        film_.U.name()
    };
}


void Foam::fv::filmCloudTransfer::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();
}


template<class Type>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::filmCloudTransfer::CloudToFilmTransferRate
(
    const Field<Type>& prop,
    const dimensionSet& dimProp
) const
{
    tmp<VolInternalField<Type>> tSu
    (
        VolInternalField<Type>::New
        (
            "Su",
            mesh(),
            dimensioned<Type>(dimProp/dimVolume/dimTime, Zero)
        )
    );

    if (prop.size())
    {
        const fvMesh& cloudMesh =
            refCast<const fvMesh>(film_.surfacePatchMap().nbrMesh());

        const label cloudPatchi =
            film_.surfacePatchMap().nbrPolyPatch().index();

        UIndirectList<Type>(tSu.ref(), film_.surfacePatch().faceCells()) =
            film_.surfacePatchMap().fromNeighbour
            (
                prop/cloudMesh.boundary()[cloudPatchi].magSf()
            );

        tSu.ref().field() /= film_.VbyA;
        tSu.ref().field() /= mesh().time().deltaTValue();
    }

    return tSu;
}


void Foam::fv::filmCloudTransfer::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    // Droplet impingement pressure
    if (fieldName == "pi")
    {
        eqn +=
            CloudToFilmTransferRate<scalar>
            (
                pressureFromCloud_,
                dimPressure*dimVolume
            );
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmCloudTransfer::addSup
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

    if (fieldName == film_.alpha.name())
    {
        eqn += CloudToFilmTransferRate<scalar>(massFromCloud_, dimMass);
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmCloudTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == film_.thermo.he().name())
    {
        eqn += CloudToFilmTransferRate<scalar>(energyFromCloud_, dimEnergy);
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmCloudTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    eqn += CloudToFilmTransferRate<vector>(momentumFromCloud_, dimMomentum);
}


void Foam::fv::filmCloudTransfer::resetFromCloudFields()
{
    const fvMesh& cloudMesh =
        refCast<const fvMesh>(film_.surfacePatchMap().nbrMesh());
    const label cloudPatchi = film_.surfacePatchMap().nbrPolyPatch().index();
    const label nCloudPatchFaces = cloudMesh.boundary()[cloudPatchi].size();

    if (massFromCloud_.size() != nCloudPatchFaces)
    {
        massFromCloud_.setSize(nCloudPatchFaces);
        momentumFromCloud_.setSize(nCloudPatchFaces);
        pressureFromCloud_.setSize(nCloudPatchFaces);
        energyFromCloud_.setSize(nCloudPatchFaces);
    }

    massFromCloud_ = 0;
    momentumFromCloud_ = Zero;
    pressureFromCloud_ = 0;
    energyFromCloud_ = 0;
}


void Foam::fv::filmCloudTransfer::parcelFromCloud
(
    const label facei,
    const scalar mass,
    const vector& momentum,
    const scalar pressure,
    const scalar energy
)
{
    massFromCloud_[facei] += mass;
    momentumFromCloud_[facei] += momentum;
    pressureFromCloud_[facei] += pressure;
    energyFromCloud_[facei] += energy;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
inline Foam::fv::filmCloudTransfer::filmToCloudTransfer
(
    const VolField<Type>& prop
) const
{
    return film_.surfacePatchMap().toNeighbour
    (
        Field<Type>(prop, film_.surfacePatch().faceCells())
    );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::filmCloudTransfer::deltaToCloud() const
{
    return filmToCloudTransfer(film_.delta);
}


void Foam::fv::filmCloudTransfer::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::filmCloudTransfer::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::filmCloudTransfer::distribute(const polyDistributionMap&)
{}


bool Foam::fv::filmCloudTransfer::movePoints()
{
    return true;
}


// ************************************************************************* //
