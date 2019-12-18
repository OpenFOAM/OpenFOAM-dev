/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::fv::constantHeatTransfer

Description
    Constant heat transfer model. htcConst [W/m^2/K] and area/volume [1/m]
    must be provided.

\*---------------------------------------------------------------------------*/

#ifndef constantHeatTransfer_H
#define constantHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class constantHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class constantHeatTransfer
:
    public interRegionHeatTransferModel
{
    // Private Data

        //- Constant heat transfer coefficient [W/m^2/K]
        autoPtr<volScalarField> htcConst_;

        //- Area per unit volume of heat exchanger [1/m]
        autoPtr<volScalarField> AoV_;


public:

    //- Runtime type information
      TypeName("constantHeatTransfer");


    // Constructors

        //- Construct from dictionary
        constantHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constantHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc();


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
