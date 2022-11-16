/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::wallBoilingModels::departureFrequencyModels::Cole

Description
    Cole correlation for bubble departure frequency.

    Reference:
    \verbatim
        Cole, R. (1960).
        A photographic study of pool boiling in the region of the critical heat
        flux.
        AIChE Journal, 6(4), 533-538.
    \endverbatim

SourceFiles
    Cole.C

\*---------------------------------------------------------------------------*/

#ifndef Cole_H
#define Cole_H

#include "departureFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureFrequencyModels
{

/*---------------------------------------------------------------------------*\
                            Class Cole Declaration
\*---------------------------------------------------------------------------*/

class Cole
:
    public departureFrequencyModel
{
public:

    //- Runtime type information
    TypeName("Cole");


    // Constructors

        //- Construct from a dictionary
        Cole(const dictionary& dict);

        //- Copy construct
        Cole(const Cole& model);

        //- Construct and return a clone
        virtual autoPtr<departureFrequencyModel> clone() const
        {
            return autoPtr<departureFrequencyModel>(new Cole(*this));
        }


    //- Destructor
    virtual ~Cole();


    // Member Functions

        //- Calculate and return the bubble departure frequency
        virtual tmp<scalarField> fDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L,
            const scalarField& dDep
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace departureFrequencyModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
