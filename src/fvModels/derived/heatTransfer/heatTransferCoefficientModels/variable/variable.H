/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::fv::heatTransferCoefficientModels::variable

Description
    Variable heat transfer model depending on local values. The Nusselt number
    is calculated as:

        \f[
            Nu = a*Re^b*Pr^c
        \f]

    And the heat transfer coefficient is calculated as:

        \f[
            htc = Nu*\kappa/L
        \f]

    Where:

    \vartable
        \kappa | Conductivity
        L      | Length scale
    \endvartable

Usage
    Example usage:
    \verbatim
    {
        heatTransferCoefficientModel variable;

        a           1;
        b           2;
        c           3;
        Pr          0.7;
        L           1e-3;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoefficientModels_variable_H
#define heatTransferCoefficientModels_variable_H

#include "heatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class variable Declaration
\*---------------------------------------------------------------------------*/

class variable
:
    public heatTransferCoefficientModel
{
    // Private Data

        //- Name of neighbour velocity field; default = U
        word UName_;

        //- Model constant
        scalar a_;

        //- Model constant
        scalar b_;

        //- Model constant
        scalar c_;

        //- Length scale
        dimensionedScalar L_;

        //- Fluid Prandtl number
        dimensionedScalar Pr_;

        //- Heat transfer coefficient [W/m^2/K]
        volScalarField htc_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("variable");


    // Constructors

        //- Construct from dictionary and mesh
        variable(const dictionary& dict, const fvMesh& mesh);

        //- Construct from dictionary and model
        variable(const dictionary& dict, const interRegionModel& model);


    //- Destructor
    virtual ~variable();


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

} // End namespace heatTransferCoefficientModels
} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
