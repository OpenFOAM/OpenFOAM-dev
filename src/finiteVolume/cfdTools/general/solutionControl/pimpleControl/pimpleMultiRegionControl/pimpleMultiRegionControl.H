/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::pimpleMultiRegionControl

Description
    PIMPLE multi-region control class. As Foam::pimpleControl, but for a multi-
    region simulation comprising PIMPLE and solid regions. More region types
    could be added.

SourceFiles
    pimpleMultiRegionControlI.H
    pimpleMultiRegionControl.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleMultiRegionControl_H
#define pimpleMultiRegionControl_H

#include "multiRegionSolutionControl.H"
#include "pimpleLoop.H"
#include "convergenceControl.H"
#include "correctorConvergenceControl.H"
#include "pimpleNoLoopControl.H"
#include "solidNoLoopControl.H"

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
protected:

    // Protected data

        //- The per-region PIMPLE controls
        PtrList<pimpleNoLoopControl> pimpleControls_;

        //- The per-region solid controls
        PtrList<solidNoLoopControl> solidControls_;


    // Protected Static Functions

        //- Get the time from the region meshes
        static const Time& time
        (
            const PtrList<fvMesh>& pimpleMeshes,
            const PtrList<fvMesh>& solidMeshes
        );


public:

    // Static data members

        //- Run-time type information
        TypeName("pimpleMultiRegionControl");


    // Constructors

        //- Construct from meshes and the algorithm name
        pimpleMultiRegionControl
        (
            PtrList<fvMesh>& pimpleMeshes,
            PtrList<fvMesh>& solidMeshes,
            const word& algorithmName="PIMPLE"
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

            //- Access the i-th solid control
            inline solidNoLoopControl& solid(const label i);

            //- Flag to indicate whether to solve the turbulence
            inline bool pimpleTurbCorr(const label i) const;

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

            //- PIMPLE loop
            bool loop();

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
