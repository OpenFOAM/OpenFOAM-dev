/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "wallDistData.H"
#include "patchDataWave.H"
#include "wallPolyPatch.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransferType>
Foam::wallDistData<TransferType>::wallDistData
(
    const Foam::fvMesh& mesh,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const bool correctWalls
)
:
    volScalarField
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimLength, great)
    ),
    cellDistFuncs(mesh),
    field_(field),
    correctWalls_(correctWalls),
    nUnset_(0)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class TransferType>
Foam::wallDistData<TransferType>::~wallDistData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransferType>
void Foam::wallDistData<TransferType>::correct()
{
    const polyMesh& mesh = cellDistFuncs::mesh();

    //
    // Fill data on wall patches with initial values
    //

    // Get patchids of walls
    labelHashSet wallPatchIDs(getPatchIDs<wallPolyPatch>());

    // Collect pointers to data on patches
    UPtrList<Field<Type>> patchData(mesh.boundaryMesh().size());

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& fieldBf = field_.boundaryFieldRef();

    forAll(fieldBf, patchi)
    {
        patchData.set(patchi, &fieldBf[patchi]);
    }

    // Do mesh wave
    patchDataWave<TransferType> wave
    (
        mesh,
        wallPatchIDs,
        patchData,
        correctWalls_
    );

    // Transfer cell values from wave into *this and field_
    transfer(wave.distance());

    field_.transfer(wave.cellData());

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& bf = boundaryFieldRef();

    // Transfer values on patches into boundaryField of *this and field_
    forAll(bf, patchi)
    {
        scalarField& waveFld = wave.patchDistance()[patchi];

        if (!isA<emptyFvPatchScalarField>(boundaryField()[patchi]))
        {
            bf[patchi].transfer(waveFld);
            Field<Type>& wavePatchData = wave.patchData()[patchi];
            fieldBf[patchi].transfer(wavePatchData);
        }
    }

    // Transfer number of unset values
    nUnset_ = wave.nUnset();
}


// ************************************************************************* //
