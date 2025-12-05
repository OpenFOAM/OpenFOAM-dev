/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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
    Foam::allGAMGProcAgglomeration

Description
    Agglomerates all processor matrices onto the master processor at the
    coarsest but one level.

    It would be more logical and useful if the agglomeration onto the master
    took place at the coarsest level, providing a means to iteratively solve the
    coarsest matrix without the need for communication, but due to a programming
    limitation of the current implementation the agglomeration takes place at
    the next finer level.  It would make more sense if this processor
    agglomerator were replaced by an \c iterativelySolveCoarsest option, similar
    to the current directSolveCoarsest option but which constructs a single
    lduMatrix which is solved using one of the CG iterative solvers.

Usage
    Example of GAMG solver settings with all processor agglomeration
    \verbatim
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;

        processorAgglomeration
        {
            agglomerator         all;
        }

        tolerance       1e-8;
        relTol          0;
    }
    \endverbatim

SourceFiles
    allGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef allGAMGProcAgglomeration_H
#define allGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class allGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class allGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        //- List of communicators
        DynamicList<label> comms_;


public:

    //- Runtime type information
    TypeName("all");


    // Constructors

        //- Construct given agglomerator and controls
        allGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        allGAMGProcAgglomeration
        (
            const allGAMGProcAgglomeration&
        ) = delete;


    //- Destructor
    virtual ~allGAMGProcAgglomeration();


    // Member Functions

        //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const allGAMGProcAgglomeration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
