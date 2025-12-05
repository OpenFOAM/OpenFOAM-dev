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
    Foam::sequentialGAMGProcAgglomeration

Description
    Sequential processor agglomeration agglomerates 2^mergeLevels sequential
    groups of processors into the minimum processor index of the group.

    Note that no account is taken of the processor connectivity so the
    processors in the agglomerated groups may or may not be connected
    depending on the decomposition ordering.

    To directly account for processor connectivity use
    Foam::pairGAMGProcAgglomeration.

Usage
    Example of GAMG solver settings with sequential processor agglomeration
    \verbatim
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;

        processorAgglomeration
        {
            agglomerator         sequential;
        }

        tolerance       1e-8;
        relTol          0;
    }
    \endverbatim

See also
    Foam::pairGAMGProcAgglomeration

SourceFiles
    sequentialGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef sequentialGAMGProcAgglomeration_H
#define sequentialGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class sequentialGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class sequentialGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        //- Agglomeration level
        const label mergeLevels_;

        DynamicList<label> comms_;


public:

    //- Runtime type information
    TypeName("sequential");


    // Constructors

        //- Construct given agglomerator and controls
        sequentialGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sequentialGAMGProcAgglomeration
        (
            const sequentialGAMGProcAgglomeration&
        ) = delete;


    //- Destructor
    virtual ~sequentialGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sequentialGAMGProcAgglomeration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
