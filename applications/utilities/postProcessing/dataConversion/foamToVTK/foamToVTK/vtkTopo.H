/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Class
    Foam::vtkTopo

Description
    Polyhedral cell decomposition for VTK.

SourceFiles
    vtkTopo.C

\*---------------------------------------------------------------------------*/

#ifndef vtkTopo_H
#define vtkTopo_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class vtkTopo Declaration
\*---------------------------------------------------------------------------*/

class vtkTopo
{
public:

   // Public Enumerations

        //- Groups of cell types retain as polyhedra
        enum class vtkPolyhedra
        {
            none,
            polyhedra,
            all
        };

        //- Names of groups of cell types retain as polyhedra
        static const NamedEnum<vtkPolyhedra, 3> vtkPolyhedraNames_;

        //- Equivalent to enumeration in "vtkCellType.h"
        enum vtkTypes
        {
            VTK_TRIANGLE   = 5,
            VTK_POLYGON    = 7,
            VTK_QUAD       = 9,
            VTK_TETRA      = 10,
            VTK_HEXAHEDRON = 12,
            VTK_WEDGE      = 13,
            VTK_PYRAMID    = 14,
            VTK_POLYHEDRON = 42
        };


private:

    // Private Data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Cell types to retain as polyhedra
        const vtkPolyhedra polyhedra_;

        //- Vertices per cell (including added cells) in vtk ordering
        labelListList vertLabels_;

        //- Cell types (including added cells) in vtk numbering
        labelList cellTypes_;

        //- ...
        labelList addPointCellLabels_;

        //- ...
        labelList superCells_;


public:

    // Constructors

        //- Construct from components
        vtkTopo(const polyMesh& mesh, const vtkPolyhedra& polyhedra);

        //- Disallow default bitwise copy construction
        vtkTopo(const vtkTopo&) = delete;


    // Member Functions

        // Access

            const labelListList& vertLabels() const
            {
                return vertLabels_;
            }

            const labelList& cellTypes() const
            {
                return cellTypes_;
            }

            const labelList& addPointCellLabels() const
            {
                return addPointCellLabels_;
            }

            const labelList& superCells() const
            {
                return superCells_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const vtkTopo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
