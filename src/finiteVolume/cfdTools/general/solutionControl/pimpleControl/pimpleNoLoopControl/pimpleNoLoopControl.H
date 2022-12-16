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
    Foam::pimpleNoLoopControl

Description
    Pimple no-loop control class. Implements various option flags, but leaves
    loop controls to the derivation or owner. Can be derived into a "full"
    pimple control or can be owned by a multi-region pimple class.

SourceFiles
    pimpleNoLoopControlI.H
    pimpleNoLoopControl.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleNoLoopControl_H
#define pimpleNoLoopControl_H

#include "pisoControl.H"
#include "singleRegionConvergenceControl.H"
#include "singleRegionCorrectorConvergenceControl.H"

#undef PISO_CONTROL
#define PIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pimpleLoop;

/*---------------------------------------------------------------------------*\
                     Class pimpleNoLoopControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleNoLoopControl
:
    public pisoControl,
    public singleRegionConvergenceControl,
    public singleRegionCorrectorConvergenceControl
{
    friend class pimpleSingleRegionControl;
    friend class pimpleMultiRegionControl;

    // Private data

        pimpleLoop* pimpleLoopPtr_;


protected:

    // Protected data

        //- Switch to move the mesh at the start of every PIMPLE outer corrector
        //  rather than the first corrector only which is the default
        bool moveMeshOuterCorrectors_;

        //- Switch to indicate whether to update the density in SIMPLE mode
        //  rather than PISO mode
        bool simpleRho_;

        //- Switch to control whether transport models are predicted
        //  just on the first pimple iteration, or on every pimple iteration
        bool transportPredictionFirst_;

        //- Switch to control whether transport models are corrected
        //  just on the final pimple iteration, or on every pimple iteration
        bool transportCorrectionFinal_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pimpleNoLoopControl");


    // Constructors

        //- Construct from a mesh, the name of the algorithm, and a reference
        //  to the pimple loop
        pimpleNoLoopControl(fvMesh& mesh, const word& algorithmName="PIMPLE");


    //- Destructor
    virtual ~pimpleNoLoopControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();

        // Access

            //- Switch to move the mesh at the start of every PIMPLE
            //  outer corrected rather than the first corrector only
            //  which is the default
            inline bool moveMeshOuterCorrectors() const;

            //- Switch to indicate whether to update the density in simple mode
            inline bool simpleRho() const;

            //- Flag to indicate whether to predict the transport models
            inline bool predictTransport() const;

            //- Flag to indicate whether to correct the transport models
            inline bool correctTransport() const;

            //- Flag to indicate the first iteration
            inline bool firstIter() const;

            //- Flag to indicate the final iteration
            inline bool finalIter() const;


        // Evolution

            //- Piso loop within outer loop
            bool correct();

            //- Non-orthogonal corrector loop
            bool correctNonOrthogonal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleNoLoopControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
