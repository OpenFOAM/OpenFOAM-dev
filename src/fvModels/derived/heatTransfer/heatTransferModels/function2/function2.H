/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fv::heatTransferModels::function2

Description
    Function2 heat transfer model. The 2D function returns the heat transfer
    coefficient as a function of the local and neighbouring velocity
    magnitudes. An area-per-unit-volume [1/m] field (AoV) must be provided in
    constant.

Usage
    Example usage:
    \verbatim
    {
        type        function2;

        AoV         1e3;

        htcFunc
        {
            type        constant;
            value       1e5;
        }
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef function2_H
#define function2_H

#include "heatTransferModel.H"
#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                         Class function2 Declaration
\*---------------------------------------------------------------------------*/

class function2
:
    public heatTransferModel
{
    // Private Data

        //- Reference to the inter region model
        const interRegionModel& model_;

        //- Name of velocity field; default = U
        word UName_;

        //- Name of neighbour velocity field; default = U
        word UNbrName_;

        //- Heat transfer coefficient function ptr
        autoPtr<Function2<scalar>> htcFunc_;

        //- Heat transfer coefficient [W/m^2/K]
        volScalarField htc_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("function2");


    // Constructors

        //- Construct from dictionary and model
        function2(const dictionary& dict, const interRegionModel& model);


    //- Destructor
    virtual ~function2();


    // Member Functions

        //- Get the heat transfer coefficient
        virtual tmp<volScalarField> htc() const
        {
            return htc_;
        }

        //- Correct the heat transfer coefficient
        virtual void correct();

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
