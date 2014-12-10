/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::surfaceInterpolation

Description
    Cell to surface interpolation scheme. Included in fvMesh.

SourceFiles
    surfaceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolation_H
#define surfaceInterpolation_H

#include "tmp.H"
#include "scalar.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class surfaceInterpolation Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolation
{
    // Private data

        // Reference to fvMesh
        const fvMesh& mesh_;

        // Demand-driven data

            //- Linear difference weighting factors
            mutable surfaceScalarField* weights_;

            //- Cell-centre difference coefficients
            mutable surfaceScalarField* deltaCoeffs_;

            //- Non-orthogonal cell-centre difference coefficients
            mutable surfaceScalarField* nonOrthDeltaCoeffs_;

            //- Non-orthogonality correction vectors
            mutable surfaceVectorField* nonOrthCorrectionVectors_;


    // Private Member Functions

        //- Construct central-differencing weighting factors
        void makeWeights() const;

        //- Construct face-gradient difference factors
        void makeDeltaCoeffs() const;

        //- Construct face-gradient difference factors
        void makeNonOrthDeltaCoeffs() const;

        //- Construct non-orthogonality correction vectors
        void makeNonOrthCorrectionVectors() const;


protected:

    // Protected Member Functions

        // Storage management

            //- Clear all geometry and addressing
            void clearOut();


public:

    // Declare name of the class and its debug switch
    ClassName("surfaceInterpolation");


    // Constructors

        //- Construct given an fvMesh
        explicit surfaceInterpolation(const fvMesh&);


    //- Destructor
    ~surfaceInterpolation();


    // Member functions

        //- Return reference to linear difference weighting factors
        const surfaceScalarField& weights() const;

        //- Return reference to cell-centre difference coefficients
        const surfaceScalarField& deltaCoeffs() const;

        //- Return reference to non-orthogonal cell-centre difference
        //  coefficients
        const surfaceScalarField& nonOrthDeltaCoeffs() const;

        //- Return reference to non-orthogonality correction vectors
        const surfaceVectorField& nonOrthCorrectionVectors() const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
