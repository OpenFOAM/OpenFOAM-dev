/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::LESModels::vanDriestDelta

Description
    Simple cube-root of cell volume delta used in incompressible LES models.

SourceFiles
    vanDriestDelta.C

\*---------------------------------------------------------------------------*/

#ifndef vanDriestDelta_H
#define vanDriestDelta_H

#include "incompressibleLESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class vanDriestDelta Declaration
\*---------------------------------------------------------------------------*/

class vanDriestDelta
:
    public LESdelta
{
    // Private data

        autoPtr<Foam::LESdelta> geometricDelta_;
        scalar kappa_;
        scalar Aplus_;
        scalar Cdelta_;
        label calcInterval_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        vanDriestDelta(const vanDriestDelta&);
        void operator=(const vanDriestDelta&);

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("vanDriest");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        vanDriestDelta(const word& name, const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~vanDriestDelta()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
