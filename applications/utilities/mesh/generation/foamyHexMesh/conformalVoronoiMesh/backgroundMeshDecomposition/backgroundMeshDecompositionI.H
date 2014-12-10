/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::backgroundMeshDecomposition::mesh() const
{
    return mesh_;
}


const Foam::indexedOctree<Foam::treeDataBPatch>&
Foam::backgroundMeshDecomposition::tree() const
{
    return bFTreePtr_();
}


const Foam::treeBoundBox&
Foam::backgroundMeshDecomposition::procBounds() const
{
    return allBackgroundMeshBounds_[Pstream::myProcNo()];
}


const Foam::labelList& Foam::backgroundMeshDecomposition::cellLevel() const
{
    return meshCutter_.cellLevel();
}


const Foam::labelList& Foam::backgroundMeshDecomposition::pointLevel() const
{
    return meshCutter_.pointLevel();
}


const Foam::decompositionMethod&
Foam::backgroundMeshDecomposition::decomposer() const
{
    return decomposerPtr_();
}


// ************************************************************************* //
