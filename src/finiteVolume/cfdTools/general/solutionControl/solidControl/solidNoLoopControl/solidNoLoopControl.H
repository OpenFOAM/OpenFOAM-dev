/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::solidNoLoopControl

Description
    Solid no-loop control class. Implements non-orthogonal and convergence
    controls, but leaves loop controls to the derivation or owner. Can be
    derived into a "full" solid control or can be owned by a multi-region class.

SourceFiles
    solidNoLoopControl.C

\*---------------------------------------------------------------------------*/

#ifndef solidNoLoopControl_H
#define solidNoLoopControl_H

#include "nonOrthogonalSolutionControl.H"
#include "pimpleLoop.H"
#include "singleRegionConvergenceControl.H"
#include "singleRegionCorrectorConvergenceControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class solidNoLoopControl Declaration
\*---------------------------------------------------------------------------*/

class solidNoLoopControl
:
    public nonOrthogonalSolutionControl,
    public singleRegionConvergenceControl,
    public singleRegionCorrectorConvergenceControl
{
protected:

    // Protected data

        //- The pimple loop
        const pimpleLoop& loop_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("solidNoLoopControl");


    // Constructors

        //- Construct from a mesh, the name of the algorithm, and a reference
        //  to the pimple loop
        solidNoLoopControl
        (
            fvMesh& mesh,
            const word& algorithmName,
            const pimpleLoop& loop
        );


    //- Destructor
    virtual ~solidNoLoopControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();

        // Access

            //- Flag to indicate whether in final state
            virtual bool isFinal() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
