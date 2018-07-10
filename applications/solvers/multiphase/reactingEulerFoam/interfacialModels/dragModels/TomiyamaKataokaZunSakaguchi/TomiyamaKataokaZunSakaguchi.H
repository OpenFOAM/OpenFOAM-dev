/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::dragModels::TomiyamaKataokaZunSakaguchi

Description
    Drag model for gas-liquid system of Tomiyama et al.

    Reference:
    \verbatim
        "Drag coefficients of single bubbles under normal and microgravity
         conditions"
        Tomiyama, A., Kataoka, I., Zun, I., Sakaguchi, T.
        JSME International Series B, Fluids and Thermal Engineering,
        Vol. 41, 1998, pp. 472-479
    \endverbatim

SourceFiles
    TomiyamaKataokaZunSakaguchi.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaKataokaZunSakaguchi_H
#define TomiyamaKataokaZunSakaguchi_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class TomiyamaKataokaZunSakaguchi Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaKataokaZunSakaguchi
:
    public dragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        //- Residual Eotvos number
        const dimensionedScalar residualEo_;


public:

    //- Runtime type information
    TypeName("TomiyamaKataokaZunSakaguchi");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaKataokaZunSakaguchi
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~TomiyamaKataokaZunSakaguchi();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
