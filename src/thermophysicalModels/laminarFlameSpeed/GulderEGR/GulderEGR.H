/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::laminarFlameSpeedModels::GulderEGR

Description
    Laminar flame speed obtained from Gulder's correlation with EGR modelling.

SourceFiles
    GulderEGR.C

\*---------------------------------------------------------------------------*/

#ifndef GulderEGR_H
#define GulderEGR_H

#include "laminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{

/*---------------------------------------------------------------------------*\
                         Class GulderEGR Declaration
\*---------------------------------------------------------------------------*/

class GulderEGR
:
    public laminarFlameSpeed
{
    // Private Data

        dictionary coeffDict_;

        scalar W_;
        scalar eta_;
        scalar xi_;
        scalar f_;
        scalar alpha_;
        scalar beta_;


    // Private Member Functions

        inline scalar SuRef
        (
            scalar phi
        ) const;

        inline scalar Su0pTphi
        (
            scalar p,
            scalar Tu,
            scalar phi,
            scalar Yres
        ) const;

        tmp<volScalarField> Su0pTphi
        (
            const volScalarField& p,
            const volScalarField& Tu,
            scalar phi
        ) const;

        tmp<volScalarField> Su0pTphi
        (
            const volScalarField& p,
            const volScalarField& Tu,
            const volScalarField& phi,
            const volScalarField& egr
        ) const;


public:

    //- Runtime type information
    TypeName("GulderEGR");


    // Constructors

        //- Construct from dictionary and psiuMulticomponentThermo
        GulderEGR
        (
            const dictionary&,
            const psiuMulticomponentThermo&
        );

        //- Disallow default bitwise copy construction
        GulderEGR(const GulderEGR&) = delete;


    //- Destructor
    virtual ~GulderEGR();


    // Member Operators

        //- Return the laminar flame speed [m/s]
        tmp<volScalarField> operator()() const;

        //- Disallow default bitwise assignment
        void operator=(const GulderEGR&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarFlameSpeedModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
