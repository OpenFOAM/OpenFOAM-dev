/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "radialActuationDiskSource.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(radialActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        radialActuationDiskSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::radialActuationDiskSource::readCoeffs()
{
    coeffs().lookup("coeffs") >> radialCoeffs_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::radialActuationDiskSource::radialActuationDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    actuationDiskSource(name, modelType, dict, mesh),
    radialCoeffs_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::radialActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const scalarField& cellsV = mesh().V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    addRadialActuationDiskAxialInertialResistance
    (
        Usource,
        set_.cells(),
        cellsV,
        geometricOneField(),
        U
    );
}


void Foam::fv::radialActuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const scalarField& cellsV = mesh().V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    addRadialActuationDiskAxialInertialResistance
    (
        Usource,
        set_.cells(),
        cellsV,
        rho,
        U
    );
}


bool Foam::fv::radialActuationDiskSource::read(const dictionary& dict)
{
    if (actuationDiskSource::read(dict))
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
