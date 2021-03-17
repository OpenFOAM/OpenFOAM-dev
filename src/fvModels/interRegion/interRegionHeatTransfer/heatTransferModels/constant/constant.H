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
    Foam::fv::heatTransferModels::constant

Description
    Constant heat transfer model. The heat transfer coefficient [W/m^2/K] (htc)
    must be provided as a value in the coefficients dictionary or as a field in
    constant.

Usage
    Example usage:
    \verbatim
    {
        type        constant;

        AoV         1e3;

        htc         1e5;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public heatTransferModel
{
    // Private Data

        //- Heat transfer coefficient [W/m^2/K]
        dimensionedScalar htc_;

        //- Heat transfer coefficient [W/m^2/K]
        autoPtr<volScalarField> htcPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary and mesh
        constant(const dictionary& dict, const fvMesh& mesh);

        //- Construct from dictionary and model
        constant(const dictionary& dict, const interRegionModel& model);


    //- Destructor
    virtual ~constant();


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
