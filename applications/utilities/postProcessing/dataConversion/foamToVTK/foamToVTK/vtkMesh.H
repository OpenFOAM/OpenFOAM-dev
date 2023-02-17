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
    Foam::vtkMesh

Description
    Encapsulation of VTK mesh data. Holds mesh or meshsubset and
    polyhedral-cell decomposition on it.

SourceFiles
    vtkMesh.C

\*---------------------------------------------------------------------------*/

#ifndef vtkMesh_H
#define vtkMesh_H

#include "vtkTopo.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;

/*---------------------------------------------------------------------------*\
                           Class vtkMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkMesh
{
    // Private Data

        //- Reference to mesh
        fvMesh& baseMesh_;

        //- Subsetting engine
        fvMeshSubset subsetter_;

        //- Cell types to retain as polyhedra
        const vtkTopo::vtkPolyhedra polyhedra_;

        //- Current cellSet (or empty)
        const word setName_;

        //- Current decomposition of topology
        mutable autoPtr<vtkTopo> topoPtr_;


public:

    // Constructors

        //- Construct from components
        vtkMesh
        (
            fvMesh& baseMesh,
            const vtkTopo::vtkPolyhedra polyhedra = vtkTopo::vtkPolyhedra::none,
            const word& setName = word::null
        );

        //- Disallow default bitwise copy construction
        vtkMesh(const vtkMesh&) = delete;


    // Member Functions

        // Access

            //- Whole mesh
            const fvMesh& baseMesh() const
            {
                return baseMesh_;
            }

            //- Subsetting engine
            const fvMeshSubset& subsetter() const
            {
                return subsetter_;
            }

            //- Check if running subMesh
            bool useSubMesh() const
            {
                return setName_.size();
            }

            //- VTK topology
            const vtkTopo& topo() const
            {
                if (topoPtr_.empty())
                {
                    topoPtr_.reset(new vtkTopo(mesh(), polyhedra_));
                }
                return topoPtr_();
            }

            //- Access either mesh or submesh
            const fvMesh& mesh() const
            {
                if (useSubMesh())
                {
                    return subsetter_.subMesh();
                }
                else
                {
                    return baseMesh_;
                }
            }

            //- Number of field cells
            label nFieldCells() const
            {
                return topo().cellTypes().size();
            }

            //- Number of field points
            label nFieldPoints() const
            {
                return mesh().nPoints() + topo().addPointCellLabels().size();
            }


        // Edit

            //- Read mesh
            polyMesh::readUpdateState readUpdate();


            //- Map volume field (does in fact do very little interpolation;
            //  just copied from fvMeshSubset)
            template<class GeoField>
            tmp<GeoField> interpolate(const GeoField& fld) const
            {
                if (useSubMesh())
                {
                    tmp<GeoField> subFld = subsetter_.interpolate(fld);
                    subFld.ref().rename(fld.name());
                    return subFld;
                }
                else
                {
                    return fld;
                }
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const vtkMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
