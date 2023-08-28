/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Foam::singleRegionCorrectorConvergenceControl

Description
    Single-region-specific derivation of the corrector convergence control
    class

SourceFiles
    singleRegionCorrectorConvergenceControl.C

\*---------------------------------------------------------------------------*/

#ifndef singleRegionCorrectorConvergenceControl_H
#define singleRegionCorrectorConvergenceControl_H

#include "correctorConvergenceControl.H"
#include "singleRegionSolutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class singleRegionCorrectorConvergenceControl Declaration
\*---------------------------------------------------------------------------*/

class singleRegionCorrectorConvergenceControl
:
    public correctorConvergenceControl
{
protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- List of residual data per field
        List<corrResidualData> corrResidualControl_;

        //- The index of the solution at the start of the corrector loop, for
        //  each field. If the field name is not in the table then the index is
        //  assumed to be zero; i.e, the first solution.
        HashTable<label> solveIndex_;


    // Protected Member Functions

        //- Read controls
        virtual bool read();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("singleRegionCorrectorConvergenceControl");


    // Constructors

        //- Construct from a solution control and the loop name
        singleRegionCorrectorConvergenceControl
        (
            const singleRegionSolutionControl& control,
            const word& loopName
        );


    //- Destructor
    virtual ~singleRegionCorrectorConvergenceControl();


    // Member Functions

        // IO

            //- Print the residual controls
            void printCorrResidualControls(const label n) const;


        // Evolution

            //- Return true if corrector residual controls are present
            virtual bool hasCorrResidualControls() const;

            //- Return true if all correction convergence checks are satisfied
            virtual bool corrCriteriaSatisfied() const;

            //- Reset the solve index in the correction residual control data
            virtual void resetCorrSolveIndex();

            //- Update the solve index in the correction residual control data
            virtual void updateCorrSolveIndex();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
