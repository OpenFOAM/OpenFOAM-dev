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

Class
    Foam::fvRegionSplit

Description
    This class enumerates disconnected regions of an fvMesh. It inherits a
    label list which it populates on construction with the region index for
    each cell in the mesh.

    Optionally additional faces can be blocked by providing a scalar field and
    value. A face is then considered blocked if the values in the adjacent
    cells (including via coupled boundaries) span the given value.

SourceFiles
    fvRegionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef fvRegionSplit_H
#define fvRegionSplit_H

#include "regionSplitBase.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class fvRegionSplit Declaration
\*---------------------------------------------------------------------------*/

class fvRegionSplit
:
    public regionSplitBase
{
    // Private Data

        //- Reference to the fvMesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Calculate region split in non-compact (global) numbering.
        void calcNonCompactRegionSplit
        (
            const globalIndex& globalFaces,
            const volScalarField& isoField,
            const scalar isoValue,
            labelList& cellRegion
        ) const;

        //- Calculate region split. Return the number of regions.
        label calcRegionSplit
        (
            const volScalarField& isoField,
            const scalar isoValue,
            labelList& cellRegion
        ) const;


public:

    //- Runtime type information
    ClassName("fvRegionSplit");


    // Constructors

        //- Construct from mesh
        fvRegionSplit(const fvMesh& mesh);

        //- Construct from mesh and iso-surface data to identify blocked faces
        fvRegionSplit
        (
            const fvMesh& mesh,
            const VolField<scalar>& isoField,
            const scalar isoValue
        );


    // Member Functions

        //- Return reference to the fvMesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
