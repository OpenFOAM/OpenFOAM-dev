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

Class
    Foam::FitData

Description
    Data for the upwinded and centred polynomial fit interpolation schemes.
    The linearCorrection_ determines whether the fit is for a corrected
    linear scheme (first two coefficients are corrections for owner and
    neighbour) or a pure upwind scheme (first coefficient is correction for
    owner; weight on face taken as 1).

SourceFiles
    FitData.C

\*---------------------------------------------------------------------------*/

#ifndef FitData_H
#define FitData_H

#include "MeshObject.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class FitData Declaration
\*---------------------------------------------------------------------------*/

template<class FitDataType, class ExtendedStencil, class Polynomial>
class FitData
:
    public MeshObject<fvMesh, MoveableMeshObject, FitDataType>
{
    // Private data

        //- The stencil the fit is based on
        const ExtendedStencil& stencil_;

        //- Is scheme correction on linear (true) or on upwind (false)
        const bool linearCorrection_;

        //- Factor the fit is allowed to deviate from the base scheme
        //  (linear or pure upwind)
        //  This limits the amount of high-order correction and increases
        //  stability on bad meshes
        const scalar linearLimitFactor_;

        //- Weights for central stencil
        const scalar centralWeight_;

        //- Dimensionality of the geometry
        const label dim_;

        //- Minimum stencil size
        const label minSize_;


protected:

        //- Find the normal direction (i) and j and k directions for face faci
        void findFaceDirs
        (
            vector& idir,        // value changed in return
            vector& jdir,        // value changed in return
            vector& kdir,        // value changed in return
            const label faci
        );

public:

    // Constructors

        //- Construct from components
        FitData
        (
            const fvMesh& mesh,
            const ExtendedStencil& stencil,
            const bool linearCorrection,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor
    virtual ~FitData()
    {}


    // Member functions

        //- Return reference to the stencil
        const ExtendedStencil& stencil() const
        {
            return stencil_;
        }

        //- Factor the fit is allowed to deviate from the base scheme
        scalar linearLimitFactor() const
        {
            return linearLimitFactor_;
        }

        //- Return weight for central stencil
        scalar centralWeight() const
        {
            return centralWeight_;
        }

        //- Dimensionality of the geometry
        label dim() const
        {
            return dim_;
        }

        //- Minimum stencil size
        label minSize() const
        {
            return minSize_;
        }

        bool linearCorrection() const
        {
            return linearCorrection_;
        }

        //- Calculate the fit for the specified face and set the coefficients
        void calcFit
        (
            scalarList& coeffsi, // coefficients to be set
            const List<point>&,  // Stencil points
            const scalar wLin,   // Weight for linear approximation (weights
                                 // nearest neighbours)
            const label faci     // Current face index
        );

        //- Calculate the fit for all the faces
        virtual void calcFit() = 0;

        //- Recalculate weights (but not stencil) when the mesh moves
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FitData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
