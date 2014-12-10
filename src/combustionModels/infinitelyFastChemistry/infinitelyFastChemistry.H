/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::combustionModels::infinitelyFastChemistry

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    infinitelyFastChemistry.C

\*---------------------------------------------------------------------------*/

#ifndef infinitelyFastChemistry_H
#define infinitelyFastChemistry_H

#include "singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class infinitelyFastChemistry Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class infinitelyFastChemistry
:
    public singleStepCombustion<CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar C_;


    // Private Member Functions

        //- Disallow copy construct
        infinitelyFastChemistry(const infinitelyFastChemistry&);

        //- Disallow default bitwise assignment
        void operator=(const infinitelyFastChemistry&);


public:

    //- Runtime type information
    TypeName("infinitelyFastChemistry");


    // Constructors

        //- Construct from components
        infinitelyFastChemistry(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~infinitelyFastChemistry();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "infinitelyFastChemistry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
