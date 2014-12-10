/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
tmp<GeometricField<Type, fvPatchField, volMesh> > cellReduce
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const CombineOp& cop
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const fvMesh& mesh = ssf.mesh();

    tmp<volFieldType> tresult
    (
        new volFieldType
        (
            IOobject
            (
                "cellReduce(" + ssf.name() + ')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", ssf.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );

    volFieldType& result = tresult();

    const labelUList& own = mesh.owner();
    const labelUList& nbr = mesh.neighbour();

    forAll(own, i)
    {
        label cellI = own[i];
        cop(result[cellI], ssf[i]);
    }
    forAll(nbr, i)
    {
        label cellI = nbr[i];
        cop(result[cellI], ssf[i]);
    }

    result.correctBoundaryConditions();

    return tresult;
}


template<class Type, class CombineOp>
tmp<GeometricField<Type, fvPatchField, volMesh> > cellReduce
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>&> tssf,
    const CombineOp& cop
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> >
        tvf(cellReduce(cop, tssf));

   tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
