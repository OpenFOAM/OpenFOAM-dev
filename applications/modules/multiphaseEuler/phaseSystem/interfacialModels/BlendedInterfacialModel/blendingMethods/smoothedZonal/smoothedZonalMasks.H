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
    Foam::blendingMethods::smoothedZonal::masks

Description
    Mesh object for storing smoothed zonal mask fields

SourceFiles
    smoothedZonalMasks.C

\*---------------------------------------------------------------------------*/

#ifndef smoothedZonalMasks_H
#define smoothedZonalMasks_H

#include "smoothedZonal.H"
#include "DemandDrivenMeshObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                       Class smoothedZonal::masks Declaration
\*---------------------------------------------------------------------------*/

class smoothedZonal::masks
:
    public DemandDrivenMeshObject<fvMesh, MoveableMeshObject, masks>
{
    // Private Data

        //- Default mask
        volScalarField defaultMask_;

        //- The indices of the zones
        const labelList zoneIndices_;

        //- Zone masks
        PtrList<volScalarField> zoneMasks_;


protected:

    friend class DemandDrivenMeshObject<fvMesh, MoveableMeshObject, masks>;


    // Protected Constructors

        //- Construct from a mesh and a list of zones
        explicit masks
        (
            const word& name,
            const fvMesh& mesh,
            const labelList& zoneIndices,
            const scalar nCellLengthScales,
            const label nIter
        );


public:

    //- Runtime type information
    TypeName("smoothedZonal::masks");


    // Constructors

        //- Disallow default bitwise copy construction
        masks(const masks&) = delete;

        //- Construct and return for the given set of zones
        static const masks& New
        (
            const fvMesh& mesh,
            const labelList& zoneIndices,
            const scalar nCellLengthScales,
            const label nIter
        );


    //- Destructor
    virtual ~masks();


    // Member Functions

        //- Update for mesh motion
        virtual bool movePoints();

        //- Access the default mask
        inline const volScalarField& operator()() const
        {
            return defaultMask_;
        }

        //- Access the zone masks
        inline const volScalarField& operator()(const label zonei) const
        {
            return zoneMasks_[zonei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
