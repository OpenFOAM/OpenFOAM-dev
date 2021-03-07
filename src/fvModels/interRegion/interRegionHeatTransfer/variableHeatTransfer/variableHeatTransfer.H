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
    Foam::fv::variableHeatTransfer

Description
    Variable heat transfer model depending on local values. The Nu number is
    calculated as:

        \f[
            Nu = a*Re^b*Pr^c
        \f]

    And the heat transfer coefficient is calculated as:

        \f[
            htc = Nu*K/ds
        \f]

    where:
    \vartable
        K   | Conductivity
        ds  | Strut diameter
    \endvartable

    An area-per-unit-volume [1/m] field (AoV) must be provided in constant.

\*---------------------------------------------------------------------------*/

#ifndef variableHeatTransfer_H
#define variableHeatTransfer_H

#include "interRegionHeatTransferModel.H"

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

        //- Model constant
        scalar a_;

        //- Model constant
        scalar b_;

        //- Model constant
        scalar c_;

        //- Strut diameter
        scalar ds_;

        //- Fluid Prandtl number
        scalar Pr_;

        //- Area per unit volume of heat exchanger
        autoPtr<volScalarField> AoV_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Correct the heat transfer coefficient
        virtual void correctHtc() const;


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


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
