/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::noneGAMGProcAgglomeration

Description
    Processor agglomeration of GAMGAgglomerations.

SourceFiles
    noneGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef noneGAMGProcAgglomeration_H
#define noneGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
              Class noneGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class noneGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noneGAMGProcAgglomeration
        (
            const noneGAMGProcAgglomeration&
        );

        //- Disallow default bitwise assignment
        void operator=(const noneGAMGProcAgglomeration&);


public:

    //- Runtime type information
    TypeName("noneGAMGProcAgglomeration");


    // Constructors

        //- Construct given agglomerator and controls
        noneGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor
    virtual ~noneGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate()
        {
            return false;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
