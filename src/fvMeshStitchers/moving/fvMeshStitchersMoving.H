/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::fvMeshStitchers::moving

Description
    Mesh stitcher for movng meshes

SourceFiles
    fvMeshStitchersMoving.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshStitchersMoving_H
#define fvMeshStitchersMoving_H

#include "fvMeshStitcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshStitchers
{

/*---------------------------------------------------------------------------*\
                Class fvMeshStitchers::moving Declaration
\*---------------------------------------------------------------------------*/

class moving
:
    public fvMeshStitcher
{
    // Private Member Functions

        //- Return whether or not this stitcher supports a changing mesh
        virtual bool changing() const
        {
            return true;
        }

        //- Correct the mesh fluxes following a conform operation
        virtual void conformCorrectMeshPhi(surfaceScalarField& phi);

        //- Initialise correction of the mesh fluxes. Create two error faces
        //  for every original face. These are initialised to have small
        //  opposing areas for stabilisation purposes. They will later be
        //  enlarged and given mesh fluxes by unconformErrorFaceCorrectMeshPhi.
        virtual void createNonConformalCorrectMeshPhiGeometry
        (
            SurfaceFieldBoundary<label>& polyFacesBf,
            surfaceVectorField& SfSf,
            surfaceVectorField& CfSf
        );

        //- Get a set of all the cells that are edge-connected to the
        //  owner-orig patches. This is the sub-set in which mesh phi gets
        //  corrected.
        labelHashSet ownerCoupledCellSet();

        //- Correct the mismatch in mesh flux across non-conformal couplings by
        //  synchronising as much as possible without affecting the total
        //  volume change. Synchronise asymmetrically, so that the owner patch
        //  takes the neighbour value. Doing this creates a volume conservation
        //  error in the owner-adjacent cells. Remove this by solving for a
        //  potential to correct the mesh fluxes across the faces that are
        //  edge-connected to the owner patches.
        //
        //  This process creates a very smooth and well behaved result. It
        //  cannot fully correct the mesh flux mismatch if the error does not
        //  globally conserve volume. In that case,
        //  unconformErrorFaceCorrectMeshPhi is also needed. Doing this
        //  correction first, however, vastly reduces the noise introduced by
        //  unconformErrorFaceCorrectMeshPhi.
        //
        void unconformInternalFaceCorrectMeshPhi(surfaceScalarField& phi);

        //- Correct the mismatch in the mesh flux by synchronising
        //  symmetrically across across non-conformal couplings. Compensate for
        //  any error that this creates by enlarging and adding mesh flux to
        //  the error faces.
        //
        //  This correction procedure can result in significant noise being
        //  introduced to the system. This is compensated by doing
        //  unconformInternalFaceCorrectMeshPhi first, which acts to smooth and
        //  reduce (and sometimes remove entirely) changes to the error faces.
        //
        void unconformErrorFaceCorrectMeshPhi
        (
            const SurfaceFieldBoundary<label>& polyFacesBf,
            surfaceVectorField& SfSf,
            surfaceVectorField& CfSf,
            surfaceScalarField& phi
        );

        //- Correct the mesh fluxes following an unconform operation
        virtual void unconformCorrectMeshPhi
        (
            const SurfaceFieldBoundary<label>& polyFacesBf,
            surfaceVectorField& SfSf,
            surfaceVectorField& CfSf,
            surfaceScalarField& phi
        );


public:

    //- Runtime type information
    TypeName("moving");


    // Constructors

        //- Construct from fvMesh
        explicit moving(fvMesh&);


    //- Destructor
    virtual ~moving();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshStitchers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
