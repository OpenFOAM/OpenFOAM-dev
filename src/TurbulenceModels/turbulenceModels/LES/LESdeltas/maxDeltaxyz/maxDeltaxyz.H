/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::LESModels::maxDeltaxyz

Description
    Delta calculated by taking the maximum distance between the cell centre
    and any face centre.  For a regular hex cell, the computed delta will
    equate to half of the cell width; accordingly, the deltaCoeff model
    coefficient should be set to 2 for this case.

SourceFiles
    maxDeltaxyz.C

\*---------------------------------------------------------------------------*/

#ifndef maxDeltaxyz_H
#define maxDeltaxyz_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class maxDeltaxyz Declaration
\*---------------------------------------------------------------------------*/

class maxDeltaxyz
:
    public LESdelta
{
    // Private data

        //- Model coefficient
        scalar deltaCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        maxDeltaxyz(const maxDeltaxyz&);
        void operator=(const maxDeltaxyz&);

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("maxDeltaxyz");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        maxDeltaxyz
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~maxDeltaxyz()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
