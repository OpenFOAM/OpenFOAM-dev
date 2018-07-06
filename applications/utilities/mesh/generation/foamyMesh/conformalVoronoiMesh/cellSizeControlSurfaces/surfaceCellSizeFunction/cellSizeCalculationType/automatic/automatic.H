/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::automatic

Description

SourceFiles
    automatic.C

\*---------------------------------------------------------------------------*/

#ifndef automatic_H
#define automatic_H

#include "cellSizeCalculationType.H"
#include "triSurfaceFields.H"
#include "PrimitivePatchInterpolation.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurfaceMesh;

/*---------------------------------------------------------------------------*\
                         Class automatic Declaration
\*---------------------------------------------------------------------------*/

class automatic
:
    public cellSizeCalculationType
{

private:

    // Private data

        //- Dictionary of coefficients for automatic cell sizing
        const dictionary& coeffsDict_;

        //- Name of the surface. Used to write the cell size field
        const fileName surfaceName_;

        const Switch readCurvature_;
        const word curvatureFile_;

        const Switch readFeatureProximity_;
        const word featureProximityFile_;

        const Switch readInternalCloseness_;
        const word internalClosenessFile_;

        //- The internalCloseness values are multiplied by the inverse
        //  of this value to get the cell size
        const scalar internalClosenessCellSizeCoeff_;

        //- The curvature values are multiplied by the inverse of this value to
        //  get the cell size
        const scalar curvatureCellSizeCoeff_;

        //- The maximum allowable sell size
        const scalar maximumCellSize_;


    // Private Member Functions

        void smoothField(triSurfaceScalarField& surf);


public:

    //- Runtime type information
    TypeName("automatic");


    // Constructors

        //- Construct from components
        automatic
        (
            const dictionary& cellSizeCalcTypeDict,
            const triSurfaceMesh& surface,
            const scalar& defaultCellSize
        );


    //- Destructor
    virtual ~automatic()
    {}


    // Member Functions

        //- Load the cell size field
        virtual tmp<triSurfacePointScalarField> load();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
