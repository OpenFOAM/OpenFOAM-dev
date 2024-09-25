/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::packingDispersionModels::Usher

Description
    Usher packing dispersion model for sewage and mineral slurries,
    both cetrifugated and filtered.

    Reference:
    \verbatim
        Usher, S. P., Studer, L. J., Wall, R. C., & Scales, P. J. (2013).
        Characterisation of dewaterability from equilibrium
        and transient centrifugation test data.
        Chemical Engineering Science, 93, 277-291.
   \endverbatim

Usage
    Example usage for mineral slurry:
    \verbatim
        packingDispersionModel Usher;

        UsherCoeffs
        {
            alphaGel        0.1;
            alphap          0.8;
            alphacp         0.8;
            alphag          0;
            alphaMax        1;
            b1              2e-3;
            sigma01         3.18663;
            n1              11;
            b2              1;
        }
    \endverbatim

    Example usage for waste water:
    \verbatim
        packingDispersionModel Usher;

        UsherCoeffs
        {
            alphaGel        0.025;
            alphap          0.09;
            alphacp         0.63;
            alphag          0.01;
            alphaMax        0.63;
            b1              0.01;
            sigma01         225.534;
            n1              6.6;
            b2              0.01;
        }
    \endverbatim

SourceFiles
    Usher.C

\*---------------------------------------------------------------------------*/

#ifndef Usher_H
#define Usher_H

#include "packingDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace packingDispersionModels
{

/*---------------------------------------------------------------------------*\
                           Class Usher Declaration
\*---------------------------------------------------------------------------*/

class Usher
:
    public packingDispersionModel
{
    // Private Data

        dimensionedScalar alphaGel_;
        dimensionedScalar alphap_;
        dimensionedScalar alphag_;
        dimensionedScalar alphacp_;
        dimensionedScalar alphaMax_;
        dimensionedScalar b1_;
        dimensionedScalar n1_;
        dimensionedScalar sigma01_;
        dimensionedScalar b2_;
        dimensionedScalar n2_;
        dimensionedScalar sigma02_;


    // Private member functions

        template<class T>
        inline T sigma1(const T& alphad) const;

        template<class T>
        inline T sigmaPrime1(const T& alphad) const;

        template<class T>
        inline T sigma2(const T& alphad) const;

        template<class T>
        inline T sigmaPrime2(const T& alphad) const;


public:

    //- Runtime type information
    TypeName("Usher");


    // Constructors

        //- Construct from components
        Usher
        (
            const dictionary& dict,
            const relativeVelocityModel& relativeVelocity
        );


    //- Destructor
    ~Usher();


    // Member Functions

        //- Return the derivative of the packing stress w.r.t. phase-fraction
        virtual tmp<volScalarField> sigmaPrime() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace packingDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
