/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

Typedefs
    CGALTriangulation3Ddefs

Description
    CGAL data structures used for 3D Delaunay meshing.

    Define CGAL_INEXACT to use Exact_predicates_inexact_constructions kernel
    otherwise the more robust but much less efficient
    Exact_predicates_exact_constructions will be used.

\*---------------------------------------------------------------------------*/

#ifndef CGALTriangulation3Ddefs_H
#define CGALTriangulation3Ddefs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CGALTriangulation3DKernel.H"

#include "indexedVertex.H"
#include "indexedCell.H"

typedef CGAL::indexedVertex<K>     Vb;
typedef CGAL::indexedCell<K>       Cb;

typedef CGAL::Compact_location     CompactLocator;
typedef CGAL::Fast_location        FastLocator;

typedef CGAL::Triangulation_data_structure_3<Vb, Cb>           Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CompactLocator> Delaunay;
typedef CGAL::Delaunay_triangulation_3<K, Tds, FastLocator>    CellSizeDelaunay;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
