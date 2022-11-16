/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    Foam::dragModels::TomiyamaAnalytic

Description
    Analytical drag model of Tomiyama et al.

    Reference:
    \verbatim
        Tomiyama, A., Kataoka, I., & Sakaguchi, T. (1995).
        Drag coefficients of bubbles (1 st Report, Drag coefficients of a single
        bubble in a stagnant liquid).
        Nihon Kikaigakkai Ronbunshu, 61, 2357-2364.
    \endverbatim

SourceFiles
    TomiyamaAnalytic.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaAnalytic_H
#define TomiyamaAnalytic_H

#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                      Class TomiyamaAnalytic Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaAnalytic
:
    public dispersedDragModel
{
    // Private Data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        //- Residual Eotvos number
        const dimensionedScalar residualEo_;

        //- Residual aspect ratio
        const dimensionedScalar residualE_;

        //- The aspect ratio model
        autoPtr<aspectRatioModel> aspectRatio_;


public:

    //- Runtime type information
    TypeName("TomiyamaAnalytic");


    // Constructors

        //- Construct from a dictionary and an interface
        TomiyamaAnalytic
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~TomiyamaAnalytic();


    // Member Functions

        // Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
