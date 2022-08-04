/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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
    Foam::pimpleMultiRegionControl

Description
    Pimple multi-region control class. As Foam::pimpleControl, but for a multi-
    region simulation comprising pimple and solid regions. More region types
    could be added.

SourceFiles
    pimpleMultiRegionControlI.H
    pimpleMultiRegionControl.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleMultiRegionControl_H
#define pimpleMultiRegionControl_H

#include "solver.H"
#include "multiRegionSolutionControl.H"
#include "pimpleLoop.H"
#include "convergenceControl.H"
#include "correctorConvergenceControl.H"
#include "pimpleNoLoopControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class pimpleMultiRegionControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleMultiRegionControl
:
    public multiRegionSolutionControl,
    public pimpleLoop,
    public convergenceControl,
    public correctorConvergenceControl
{
    // Private data

        //- The per-region pimple controls
        UPtrList<pimpleNoLoopControl> pimpleControls_;

        //- Maximum number of energy correctors
        label nEcorr_;

        //- Current energy corrector
        label Ecorr_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pimpleMultiRegionControl");


    // Constructors

        //- Construct from time and the list of solvers
        pimpleMultiRegionControl
        (
            const Time& runTime,
            PtrList<solver>& solvers,
            const word& algorithmName = "PIMPLE"
        );


    //- Destructor
    virtual ~pimpleMultiRegionControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();

        // Access

            //- Access the i-th pimple control
            inline pimpleNoLoopControl& pimple(const label i);

        // Evolution

            //- Return true if residual controls are present
            virtual bool hasResidualControls() const;

            //- Return true if corrector residual controls are present
            virtual bool hasCorrResidualControls() const;

            //- Return true if all convergence checks are satisfied
            virtual bool criteriaSatisfied() const;

            //- Return true if all correction convergence checks are satisfied
            virtual bool corrCriteriaSatisfied() const;

            //- Reset the solve index in the correction residual control data
            virtual void resetCorrSolveIndex();

            //- Update the solve index in the correction residual control data
            virtual void updateCorrSolveIndex();

            //- Pimple loop
            bool loop();

            //- Energy corrector loop
            bool correctEnergy();

            //- Time run loop
            bool run(Time& time);

            //- Time loop loop
            bool loop(Time& time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleMultiRegionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
