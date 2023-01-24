/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::wallBoilingModels::partitioningModels::phaseFraction

Description
    Wall heat-flux partitioned according to the phase volume fraction.

SourceFiles
    phaseFraction.C

\*---------------------------------------------------------------------------*/

#ifndef phaseFraction_H
#define phaseFraction_H

#include "partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                        Class phaseFraction Declaration
\*---------------------------------------------------------------------------*/

class phaseFraction
:
    public partitioningModel
{

public:

    //- Runtime type information
    TypeName("phaseFraction");


    // Constructors

        //- Construct from a dictionary
        phaseFraction(const dictionary& dict);

        //- Copy construct
        phaseFraction(const phaseFraction& model);

        //- Construct and return a clone
        virtual autoPtr<partitioningModel> clone() const
        {
            return autoPtr<partitioningModel>(new phaseFraction(*this));
        }


    //- Destructor
    virtual ~phaseFraction();


    // Member Functions

        //- Return the wet fraction
        virtual tmp<scalarField> wetFraction
        (
            const scalarField& alphaLiquid
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
