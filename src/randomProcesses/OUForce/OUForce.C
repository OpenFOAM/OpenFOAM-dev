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

#include "OUForce.H"
#include "fvMatrices.H"
#include "fft.H"
#include "writeEk.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(OUForce, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        OUForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::OUForce::readCoeffs(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("U", "U");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::OUForce::OUForce
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    UName_(word::null),
    K_(mesh),
    forceGen_(K_, mesh.time().deltaTValue(), dict)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::OUForce::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::OUForce::addSup
(
    const volVectorField& field,
    fvMatrix<vector>& eqn
) const
{
    eqn += volVectorField::Internal
    (
        IOobject
        (
            "OUForce",
            mesh().time().name(),
            mesh()
        ),
        mesh(),
        dimAcceleration,
        ReImSum
        (
            fft::reverseTransform
            (
                (K_/(mag(K_) + 1.0e-6))
              ^ forceGen_.newField(mesh().time().deltaTValue()), K_.nn()
            )
        )
    );
}


bool Foam::fv::OUForce::movePoints()
{
    return true;
}


void Foam::fv::OUForce::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::OUForce::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::OUForce::distribute(const polyDistributionMap&)
{}


bool Foam::fv::OUForce::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::fv::OUForce::write(const bool write) const
{
    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);
    writeEk(U, K_);

    return true;
}


// ************************************************************************* //
