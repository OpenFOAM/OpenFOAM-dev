/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "fviCellReduce.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
Foam::tmp<Foam::VolInternalField<Type>>
Foam::fvi::cellReduce
(
    const SurfaceField<Type>& ssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    const fvMesh& mesh = ssf.mesh()();

    tmp<VolInternalField<Type>> tresult
    (
        VolInternalField<Type>::New
        (
            "cellReduce(" + ssf.name() + ')',
            mesh,
            dimensioned<Type>
            (
                ssf.dimensions(),
                nullValue
            )
        )
    );

    VolInternalField<Type>& result = tresult.ref();

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

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::VolInternalField<Type>>
Foam::fvi::cellReduce
(
    const tmp<SurfaceField<Type>>& tssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    tmp<VolInternalField<Type>> tvf
    (
        cellReduce
        (
            tssf(),
            cop,
            nullValue
        )
    );

    tssf.clear();

    return tvf;
}


// ************************************************************************* //
