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
    Foam::XiProfiles::cubic

Description
    Simple cubic profile for the equilibrium flame wrinkling \c XiEq
    to give a plausible distribution.

    References:
    \verbatim
        Weller, H. G. (1993).
        The development of a new flame area combustion model
        using conditional averaging.
        Thermo-fluids section report TF 9307.

        Weller, H. G., Tabor, G., Gosman, A. D., & Fureby, C. (1998, January).
        Application of a flame-wrinkling LES combustion model
        to a turbulent mixing layer.
        In Symposium (International) on combustion
        (Vol. 27, No. 1, pp. 899-907). Elsevier.
    \endverbatim

SourceFiles
    cubic.C

\*---------------------------------------------------------------------------*/

#ifndef cubic_H
#define cubic_H

#include "XiProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiProfiles
{

/*---------------------------------------------------------------------------*\
                           Class cubic Declaration
\*---------------------------------------------------------------------------*/

class cubic
:
    public XiProfile
{
    // Private Data

        dimensionedScalar XiShapeCoeff_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("cubic");


    // Constructors

        //- Construct from components
        cubic
        (
            const dictionary& dict,
            const volScalarField& b
        );

        //- Disallow default bitwise copy construction
        cubic(const cubic&) = delete;


    //- Destructor
    virtual ~cubic();


    // Member Functions

        //- Return the flame-wrinkling profile
        virtual tmp<volScalarField> profile() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cubic&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiProfiles
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
