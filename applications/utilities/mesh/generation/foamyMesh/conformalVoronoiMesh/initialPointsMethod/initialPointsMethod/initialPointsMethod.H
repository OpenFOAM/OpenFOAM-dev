/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::initialPointsMethod

Description
    Abstract base class for generating initial points for a conformalVoronoiMesh

SourceFiles
    initialPointsMethod.C

\*---------------------------------------------------------------------------*/

#ifndef initialPointsMethod_H
#define initialPointsMethod_H

#include "point.H"
#include "conformalVoronoiMesh.H"
#include "backgroundMeshDecomposition.H"
#include "dictionary.H"
#include "Random.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class initialPointsMethod Declaration
\*---------------------------------------------------------------------------*/

class initialPointsMethod
:
    public dictionary
{

protected:

    // Protected data

        const Time& runTime_;

        Random& rndGen_;

        const conformationSurfaces& geometryToConformTo_;

        const cellShapeControl& cellShapeControls_;

        const autoPtr<backgroundMeshDecomposition>& decomposition_;

        //- Method details dictionary
        dictionary detailsDict_;

        //- Only allow the placement of initial points that are within the
        //  surfaces to be meshed by minimumSurfaceDistanceCoeff multiplied by
        //  the local target cell size.  Store square of value.
        scalar minimumSurfaceDistanceCoeffSqr_;

        Switch fixInitialPoints_;


public:

    //- Runtime type information
    TypeName("initialPointsMethod");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            initialPointsMethod,
            dictionary,
            (
                const dictionary& initialPointsDict,
                const Time& runTime,
                Random& rndGen,
                const conformationSurfaces& geometryToConformTo,
                const cellShapeControl& cellShapeControls,
                const autoPtr<backgroundMeshDecomposition>& decomposition
            ),
            (
                initialPointsDict,
                runTime,
                rndGen,
                geometryToConformTo,
                cellShapeControls,
                decomposition
            )
        );


    // Constructors

        //- Construct from components
        initialPointsMethod
        (
            const word& type,
            const dictionary& initialPointsDict,
            const Time& runTime,
            Random& rndGen,
            const conformationSurfaces& geometryToConformTo,
            const cellShapeControl& cellShapeControls,
            const autoPtr<backgroundMeshDecomposition>& decomposition
        );

        //- Disallow default bitwise copy construction
        initialPointsMethod(const initialPointsMethod&) = delete;


    // Selectors

        //- Return a reference to the selected initialPointsMethod
        static autoPtr<initialPointsMethod> New
        (
            const dictionary& initialPointsDict,
            const Time& runTime,
            Random& rndGen,
            const conformationSurfaces& geometryToConformTo,
            const cellShapeControl& cellShapeControls,
            const autoPtr<backgroundMeshDecomposition>& decomposition
        );


    //- Destructor
    virtual ~initialPointsMethod();


    // Member Functions

        // Access

            const Time& time() const
            {
                return runTime_;
            }

            Random& rndGen() const
            {
                return rndGen_;
            }

            const conformationSurfaces& geometryToConformTo() const
            {
                return geometryToConformTo_;
            }

            const cellShapeControl& cellShapeControls() const
            {
                return cellShapeControls_;
            }

            const backgroundMeshDecomposition& decomposition() const
            {
                return decomposition_;
            }

            //- Const access to the details dictionary
            const dictionary& detailsDict() const
            {
                return detailsDict_;
            }

            Switch fixInitialPoints() const
            {
                return fixInitialPoints_;
            }


        // Queries

            //- Return the initial points for the conformalVoronoiMesh
            virtual List<Vb::Point> initialPoints() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const initialPointsMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
