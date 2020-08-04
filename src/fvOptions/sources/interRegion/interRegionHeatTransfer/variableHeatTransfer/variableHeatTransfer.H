/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::fv::variableHeatTransfer

Description
    Variable heat transfer model depending on local values. The area of contact
    between regions (area) must be provided. The Nu number is calculated as:

        Nu = a*pow(Re, b)*pow(Pr, c)

    and the heat transfer coefficient as:

        htc = Nu*K/ds

    where:
        K is the heat conduction
        ds is the strut diameter

\*---------------------------------------------------------------------------*/

#ifndef variableHeatTransfer_H
#define variableHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class variableHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class variableHeatTransfer
:
    public interRegionHeatTransferModel
{
    // Private Data

        //- Name of neighbour velocity field; default = U
        word UNbrName_;

        //- Model constants
        scalar a_;
        scalar b_;
        scalar c_;

        //- Strut diameter
        scalar ds_;

        //- Fluid Prandtl number
        scalar Pr_;

        //- Area per unit volume of heat exchanger
        autoPtr<volScalarField> AoV_;


public:

    //- Runtime type information
    TypeName("variableHeatTransfer");


    // Constructors

        //- Construct from dictionary
        variableHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~variableHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc() const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
