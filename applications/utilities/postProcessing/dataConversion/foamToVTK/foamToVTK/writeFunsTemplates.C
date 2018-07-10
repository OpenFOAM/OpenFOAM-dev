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

#include "writeFuns.H"
#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Store List in dest
template<class Type>
void Foam::writeFuns::insert
(
    const List<Type>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(source, i)
    {
        insert(source[i], dest);
    }
}


//// Store List (indexed through map) in dest
//template<class Type>
//void Foam::writeFuns::insert
//(
//    const labelList& map,
//    const List<Type>& source,
//    DynamicList<floatScalar>& dest
//)
//{
//    forAll(map, i)
//    {
//        insert(source[map[i]], dest);
//    }
//}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const GeometricField<Type, fvPatchField, volMesh>& vvf,
    const vtkMesh& vMesh
)
{
    const fvMesh& mesh = vMesh.mesh();

    const labelList& superCells = vMesh.topo().superCells();

    label nValues = mesh.nCells() + superCells.size();

    os  << vvf.name() << ' ' << pTraits<Type>::nComponents << ' '
        << nValues << " float" << std::endl;

    DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nValues);

    insert(vvf.primitiveField(), fField);

    forAll(superCells, superCelli)
    {
        label origCelli = superCells[superCelli];

        insert(vvf[origCelli], fField);
    }
    write(os, binary, fField);
}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const GeometricField<Type, pointPatchField, pointMesh>& pvf,
    const vtkMesh& vMesh
)
{
    const fvMesh& mesh = vMesh.mesh();
    const vtkTopo& topo = vMesh.topo();

    const labelList& addPointCellLabels = topo.addPointCellLabels();
    const label nTotPoints = mesh.nPoints() + addPointCellLabels.size();

    os  << pvf.name() << ' ' << pTraits<Type>::nComponents << ' '
        << nTotPoints << " float" << std::endl;

    DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nTotPoints);

    insert(pvf, fField);

    forAll(addPointCellLabels, api)
    {
        label origCelli = addPointCellLabels[api];

        insert(interpolatePointToCell(pvf, origCelli), fField);
    }
    write(os, binary, fField);
}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const GeometricField<Type, fvPatchField, volMesh>& vvf,
    const GeometricField<Type, pointPatchField, pointMesh>& pvf,
    const vtkMesh& vMesh
)
{
    const fvMesh& mesh = vMesh.mesh();
    const vtkTopo& topo = vMesh.topo();

    const labelList& addPointCellLabels = topo.addPointCellLabels();
    const label nTotPoints = mesh.nPoints() + addPointCellLabels.size();

    os  << vvf.name() << ' ' << pTraits<Type>::nComponents << ' '
        << nTotPoints << " float" << std::endl;

    DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nTotPoints);

    insert(pvf, fField);

    forAll(addPointCellLabels, api)
    {
        label origCelli = addPointCellLabels[api];

        insert(vvf[origCelli], fField);
    }
    write(os, binary, fField);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const PtrList<GeometricField<Type, PatchField, GeoMesh>>& flds,
    const vtkMesh& vMesh
)
{
    forAll(flds, i)
    {
        write(os, binary, flds[i], vMesh);
    }
}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const volPointInterpolation& pInterp,
    const PtrList<GeometricField<Type, fvPatchField, volMesh>>& flds,
    const vtkMesh& vMesh
)
{
    forAll(flds, i)
    {
        write(os, binary, flds[i], pInterp.interpolate(flds[i])(), vMesh);
    }
}


// ************************************************************************* //
