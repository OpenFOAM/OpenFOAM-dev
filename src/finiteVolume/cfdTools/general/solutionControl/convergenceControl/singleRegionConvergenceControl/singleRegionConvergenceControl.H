/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::singleRegionConvergenceControl

Description
    Single-region-specific derivation of the convergence control class

SourceFiles
    singleRegionConvergenceControl.C

\*---------------------------------------------------------------------------*/

#ifndef singleRegionConvergenceControl_H
#define singleRegionConvergenceControl_H

#include "convergenceControl.H"
#include "singleRegionSolutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class singleRegionConvergenceControl Declaration
\*---------------------------------------------------------------------------*/

class singleRegionConvergenceControl
:
    public convergenceControl
{
protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- List of residual data per field
        List<residualData> residualControl_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("singleRegionConvergenceControl");


    // Constructors

        //- Construct from a solution control
        singleRegionConvergenceControl
        (
            const singleRegionSolutionControl& control
        );


    //- Destructor
    virtual ~singleRegionConvergenceControl();


    // Member Functions

        // IO

            //- Read residual controls
            bool readResidualControls();

            //- Print the residual controls
            void printResidualControls() const;


        // Evolution

            //- Return true if residual controls are present
            virtual bool hasResidualControls() const;

            //- Return true if all convergence checks are satisfied
            virtual bool criteriaSatisfied() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
