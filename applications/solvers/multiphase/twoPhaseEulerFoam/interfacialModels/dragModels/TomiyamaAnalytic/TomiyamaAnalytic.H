/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Analytic drag model of Tomiyama et al.

    Reference:
    \verbatim
        "Drag Coefficients of Bubbles. 1st Report. Drag Coefficients of a
        Single Bubble in a Stagnant Liquid."
        A Tomiyama, I Kataoka, and T Sakaguchi
        Nippon Kikai Gakkai Ronbunshu
        Volume 61, Issue 587, Pages 2357-2364, 1995
    \endverbatim

SourceFiles
    TomiyamaAnalytic.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaAnalytic_H
#define TomiyamaAnalytic_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class TomiyamaAnalytic Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaAnalytic
:
    public dragModel
{
private:

    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        //- Residual Eotvos number
        const dimensionedScalar residualEo_;

        //- Residual aspect ratio
        const dimensionedScalar residualE_;


public:

    //- Runtime type information
    TypeName("TomiyamaAnalytic");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaAnalytic
        (
            const dictionary& dict,
            const phasePair& pair,
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
