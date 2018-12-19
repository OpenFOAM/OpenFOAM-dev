/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "fvcCellReduce.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvc::cellReduce
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const fvMesh& mesh = ssf.mesh();

    tmp<volFieldType> tresult
    (
        volFieldType::New
        (
            "cellReduce(" + ssf.name() + ')',
            mesh,
            dimensioned<Type>
            (
                ssf.dimensions(),
                nullValue
            ),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );

    volFieldType& result = tresult.ref();

    const labelUList& own = mesh.owner();
    const labelUList& nbr = mesh.neighbour();

    forAll(own, i)
    {
        label celli = own[i];
        cop(result[celli], ssf[i]);
    }
    forAll(nbr, i)
    {
        label celli = nbr[i];
        cop(result[celli], ssf[i]);
    }

    result.correctBoundaryConditions();

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvc::cellReduce
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        cellReduce
        (
            tssf,
            cop,
            nullValue
        )
    );

    tssf.clear();

    return tvf;
}


// ************************************************************************* //
