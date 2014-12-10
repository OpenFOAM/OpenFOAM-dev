/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::PrandtlDelta

Description
    Simple cube-root of cell volume delta used in LES models.

SourceFiles
    PrandtlDelta.C

\*---------------------------------------------------------------------------*/

#ifndef PrandtlDelta_H
#define PrandtlDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PrandtlDelta Declaration
\*---------------------------------------------------------------------------*/

class PrandtlDelta
:
    public LESdelta
{
    // Private data

        autoPtr<LESdelta> geometricDelta_;
        scalar kappa_;
        scalar Cdelta_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        PrandtlDelta(const PrandtlDelta&);
        void operator=(const PrandtlDelta&);

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("Prandtl");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        PrandtlDelta
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );


    //- Destructor
    virtual ~PrandtlDelta()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
