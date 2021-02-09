/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::curvatureSeparation

Description
    Curvature film separation model

    Assesses film curvature via the mesh geometry and calculates a force
    balance of the form:

        F_sum = F_inertial + F_body + F_surface

    If F_sum < 0, the film separates. Similarly, if F_sum > 0 the film will
    remain attached.

    Based on description given by
        Owen and D. J. Ryley. The flow of thin liquid films around corners.
        International Journal of Multiphase Flow, 11(1):51-62, 1985.


SourceFiles
    curvatureSeparation.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureSeparation_H
#define curvatureSeparation_H

#include "ejectionModel.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class curvatureSeparation Declaration
\*---------------------------------------------------------------------------*/

class curvatureSeparation
:
    public ejectionModel
{
protected:

    // Protected data

        //- Gradient of surface normals
        volTensorField gradNHat_;

        //- Minimum gravity driven film thickness (non-dimensionalised delta/R1)
        scalar deltaByR1Min_;

        //- Stable film thickness - drips only formed if thickness
        //  exceeds this threshold value
        scalar deltaStable_;

        //- List of radii for patches - if patch not defined, radius
        // calculated based on mesh geometry
        List<Tuple2<label, scalar>> definedPatchRadii_;

        //- Magnitude of gravity vector
        scalar magG_;

        //- Direction of gravity vector
        vector gHat_;


    // Protected Member Functions

        //- Calculate local (inverse) radius of curvature
        tmp<volScalarField> calcInvR1(const volVectorField& U) const;

        //- Calculate the cosine of the angle between gravity vector and
        //  cell out flow direction
        tmp<scalarField> calcCosAngle(const surfaceScalarField& phi) const;


public:

    //- Runtime type information
    TypeName("curvatureSeparation");


    // Constructors

        //- Construct from surface film model
        curvatureSeparation
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        curvatureSeparation(const curvatureSeparation&) = delete;


    //- Destructor
    virtual ~curvatureSeparation();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                scalarField& massToEject,
                scalarField& diameterToEject
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const curvatureSeparation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
