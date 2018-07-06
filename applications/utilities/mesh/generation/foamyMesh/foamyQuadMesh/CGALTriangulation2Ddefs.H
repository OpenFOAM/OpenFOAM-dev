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

Typedefs
    CGALTriangulation2Ddefs

Description
    CGAL data structures used for 2D Delaunay meshing.

    Define CGAL_INEXACT to use Exact_predicates_inexact_constructions kernel
    otherwise the more robust but much less efficient
    Exact_predicates_exact_constructions will be used.

    Define CGAL_HIERARCHY to use hierarchical Delaunay triangulation which is
    faster but uses more memory than the standard Delaunay triangulation.

\*---------------------------------------------------------------------------*/

#ifndef CGALTriangulation2Ddefs_H
#define CGALTriangulation2Ddefs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CGALTriangulation2DKernel.H"

#include "indexedVertex.H"
#include "indexedFace.H"

typedef CGAL::indexedVertex<K>                           Vb;
typedef CGAL::indexedFace<K>                             Fb;

#ifdef CGAL_HIERARCHY

    // Data structures for hierarchical Delaunay triangulation which is more
    // efficient but also uses more storage
    #include "CGAL/Triangulation_hierarchy_2.h"

    typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vb>  Vbh;
    typedef CGAL::Triangulation_data_structure_2<Vbh, Fb>    Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>           DT;
    typedef CGAL::Triangulation_hierarchy_2<DT>              Delaunay;

#else

    // Data structures for standard Delaunay triangulation
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>     Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>           Delaunay;

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
