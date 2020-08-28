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
    Foam::masterCoarsestGAMGProcAgglomeration

Description
    Processor agglomeration of GAMGAgglomerations.

SourceFiles
    masterCoarsestGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef masterCoarsestGAMGProcAgglomeration_H
#define masterCoarsestGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
             Class masterCoarsestGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class masterCoarsestGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        DynamicList<label> comms_;


public:

    //- Runtime type information
    TypeName("masterCoarsest");


    // Constructors

        //- Construct given agglomerator and controls
        masterCoarsestGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );

        //- Disallow default bitwise copy construction
        masterCoarsestGAMGProcAgglomeration
        (
            const masterCoarsestGAMGProcAgglomeration&
        );


    //- Destructor
    virtual ~masterCoarsestGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();

        //- Disallow default bitwise assignment


    // Member Operators

        void operator=(const masterCoarsestGAMGProcAgglomeration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
