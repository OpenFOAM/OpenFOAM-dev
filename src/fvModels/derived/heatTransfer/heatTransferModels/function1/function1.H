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
    Foam::fv::heatTransferModels::function1

Description
    Function1 heat transfer model. The 1D function returns the heat transfer
    coefficient as a function of the local velocity magnitude.

Usage
    Example usage:
    \verbatim
    {
        type        function1;

        AoV         1e3;

        htcFunc
        {
            type        constant;
            value       1e5;
        }
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef function1_H
#define function1_H

#include "heatTransferModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                           Class function1 Declaration
\*---------------------------------------------------------------------------*/

class function1
:
    public heatTransferModel
{
    // Private Data

        //- Name of velocity field; default = U
        word UName_;

        //- Heat transfer coefficient [W/m^2/K]
        autoPtr<Function1<scalar>> htcFunc_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("function1");


    // Constructors

        //- Construct from dictionary and mesh
        function1(const dictionary& dict, const fvMesh& mesh);

        //- Construct from dictionary and model
        function1(const dictionary& dict, const interRegionModel& model);


    //- Destructor
    virtual ~function1();


    // Member Functions

        //- Get the heat transfer coefficient
        virtual tmp<volScalarField> htc() const;

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
