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
    Foam::fv::tabulatedHeatTransfer

Description
    Tabulated heat transfer model. The heat exchange area per unit volume
    must be provided.  The 2D table returns the heat transfer coefficient
    by querying the local and neighbour region velocities

\*---------------------------------------------------------------------------*/

#ifndef tabulatedHeatTransfer_H
#define tabulatedHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class tabulatedHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedHeatTransfer
:
    public interRegionHeatTransferModel
{
    // Private Data

        //- Name of velocity field; default = U
        word UName_;

        //- Name of neighbour velocity field; default = U
        word UNbrName_;

        //- 2D look up table
        mutable autoPtr<interpolation2DTable<scalar>> hTable_;

        //- Area per unit volume of heat exchanger
        mutable autoPtr<volScalarField> AoV_;

        //- Heat transfer coefficient table
        const interpolation2DTable<scalar>& hTable() const;

        //- Field of area divided by volume
        const volScalarField& AoV() const;

        //- Start time name
        const word startTimeName_;


public:

    //- Runtime type information
        TypeName("tabulatedHeatTransfer");


    // Constructors

        //- Construct from dictionary
        tabulatedHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~tabulatedHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc() const;


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
