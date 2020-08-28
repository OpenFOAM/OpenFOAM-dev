/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::eagerGAMGProcAgglomeration

Description
    'Eager' processor agglomeration of GAMGAgglomerations: at every
    level agglomerates 'mergeLevels' number of processors onto the
    minimum processor number.

SourceFiles
    eagerGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef eagerGAMGProcAgglomeration_H
#define eagerGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
                 Class eagerGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class eagerGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        //- Agglpmeration level
        const label mergeLevels_;

        DynamicList<label> comms_;


public:

    //- Runtime type information
    TypeName("eager");


    // Constructors

        //- Construct given agglomerator and controls
        eagerGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );

        //- Disallow default bitwise copy construction
        eagerGAMGProcAgglomeration
        (
            const eagerGAMGProcAgglomeration&
        );


    //- Destructor
    virtual ~eagerGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const eagerGAMGProcAgglomeration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
