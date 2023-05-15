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
    Foam::filmEjectionModels::curvatureSeparation

Description
    Curvature induced separation film to cloud ejection transfer model

    Assesses film curvature via the mesh geometry and calculates a force
    balance of the form:

        F_sum = F_inertial + F_body + F_surface_tension

    If F_sum < 0, the film separates and is transferred to the cloud
    if F_sum >= 0 the film remains attached.

    Reference:
    \verbatim
        Owen, I., & Ryley, D. J. (1985).
        The flow of thin liquid films around corners.
        International journal of multiphase flow, 11(1), 51-62.
    \endverbatim

Usage
    Example usage:
    \verbatim
    filmCloudTransfer
    {
        type    filmCloudTransfer;

        libs    ("libfilmCloudTransfer.so");

        ejection
        {
            model   curvatureSeparation;

            deltaStable 5e-4;
        }
    }
    \endverbatim

SourceFiles
    curvatureSeparation.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureSeparation_H
#define curvatureSeparation_H

#include "ejectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmEjectionModels
{

/*---------------------------------------------------------------------------*\
                     Class curvatureSeparation Declaration
\*---------------------------------------------------------------------------*/

class curvatureSeparation
:
    public ejectionModel
{
    // Private Member Data

        //- Gradient of film wall normals
        tensorField gradNHat_;

        //- Minimum gravity driven film thickness (non-dimensionalised delta/R1)
        scalar deltaByR1Min_;

        //- Stable film thickness - drips only formed if thickness
        //  exceeds this threshold value
        scalar deltaStable_;

        //- List of radii for patches
        //  if patch not listed the radius is calculated from the face normals
        List<Tuple2<label, scalar>> patchRadii_;


    // Protected Member Functions

        //- Calculate local (inverse) radius of curvature
        tmp<scalarField> calcInvR1(const volVectorField& U) const;

        //- Calculate the cosine of the angle between gravity vector and
        //  cell out flow direction
        tmp<scalarField> calcCosAngle(const surfaceScalarField& phi) const;


public:

    //- Runtime type information
    TypeName("curvatureSeparation");


    // Constructors

        //- Construct from dictionary and film model
        curvatureSeparation
        (
            const dictionary& dict,
            const solvers::isothermalFilm& film
        );

        //- Disallow default bitwise copy construction
        curvatureSeparation(const curvatureSeparation&) = delete;


    //- Destructor
    virtual ~curvatureSeparation();


    // Member Functions

        //- Correct
        virtual void correct();


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const curvatureSeparation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmEjectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
