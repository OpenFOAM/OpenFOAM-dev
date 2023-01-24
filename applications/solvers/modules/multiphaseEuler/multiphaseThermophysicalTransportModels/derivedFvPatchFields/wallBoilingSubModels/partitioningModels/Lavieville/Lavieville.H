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
    Foam::wallBoilingModels::partitioningModels::Lavieville

Description
    Lavieville wall heat flux partitioning model.

    Model parameters:
        alphaCrit: critical liquid fraction

    Reference:
    \verbatim
        Lavieville, J., Quemerais, E., Mimouni, S., Boucker, M., &
        Mechitoua, N. (2006).
        NEPTUNE CFD V1. 0 theory manual.
        NEPTUNE report Nept_2004_L1, 2(3).
    \endverbatim

SourceFiles
    Lavieville.C

\*---------------------------------------------------------------------------*/

#ifndef Lavieville_H
#define Lavieville_H

#include "partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                         Class Lavieville Declaration
\*---------------------------------------------------------------------------*/

class Lavieville
:
    public partitioningModel
{
    // Private Data

        //- Critical liquid fraction
        scalar alphaCrit_;


public:

    //- Runtime type information
    TypeName("Lavieville");


    // Constructors

        //- Construct from a dictionary
        Lavieville(const dictionary& dict);

        //- Copy construct
        Lavieville(const Lavieville& model);

        //- Construct and return a clone
        virtual autoPtr<partitioningModel> clone() const
        {
            return autoPtr<partitioningModel>(new Lavieville(*this));
        }


    //- Destructor
    virtual ~Lavieville();


    // Member Functions

        //- Return the wet fraction
        virtual tmp<scalarField> wetFraction
        (
            const scalarField& alphaLiquid
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
