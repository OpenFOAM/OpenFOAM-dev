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
    Foam::fvMeshMovers::interpolator

Description
    Interpolates pre-specified motion specified as a set of pointVectorFields.

    The motion can be provided either as a set of displacement or position
    fields and the entry \c displacement specified accordingly.

Usage
    Example:
    \verbatim
    mover   interpolator;

    displacementLaplacianCoeffs
    {
        field               wantedDisplacement;
        displacement        yes;
        interpolationScheme linear;
    }
    \endverbatim

    This will scan the case for \c wantedDisplacement \c pointVectorFields in
    the time directories and interpolate those in time (using \c linear
    interpolation) to obtain the current displacement.  The advantage of
    specifying displacement in this way is that it automatically works in
    parallel using \c decomposePar to decompose the set of \c pointVectorFields
    provided.

SourceFiles
    interpolator.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMoversInterpolator_H
#define fvMeshMoversInterpolator_H

#include "fvMeshMover.H"
#include "dynamicMeshPointInterpolator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                  Class interpolator Declaration
\*---------------------------------------------------------------------------*/

class interpolator
:
    public fvMeshMover
{
    // Private Data

        dictionary meshCoeffs_;

        dynamicMeshPointInterpolator pointInterpolator_;

        const Switch displacement_;

        //- Starting points
        autoPtr<pointVectorField> points0_;

        //- Optional list of vectorFields to update for mesh motion
        //  For modern solvers using Uf and correctPhi to update the flux
        //  after motion it is not necessary to specify a "velocityFields" list
        velocityMotionCorrection velocityMotionCorrection_;


public:

    //- Runtime type information
    TypeName("interpolator");


    // Constructors

        //- Construct from fvMesh
        interpolator(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        interpolator(const interpolator&) = delete;


    //- Destructor
    ~interpolator();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interpolator&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
