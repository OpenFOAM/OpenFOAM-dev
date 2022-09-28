/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::compressible::cavitationModels::Merkle

Description
    Merkle cavitation model.

    Reference:
    \verbatim
        C. L. Merkle, J. Feng, and P. E. O. Buelow,
        "Computational modeling of the dynamics of sheet cavitation",
        in Proceedings Third International Symposium on Cavitation
        Grenoble, France 1998.
   \endverbatim

SourceFiles
    Merkle.C

\*---------------------------------------------------------------------------*/

#ifndef Merkle_H
#define Merkle_H

#include "cavitationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace cavitationModels
{

/*---------------------------------------------------------------------------*\
                                 Class Merkle
\*---------------------------------------------------------------------------*/

class Merkle
:
    public cavitationModel
{
    // Private Data

        dimensionedScalar UInf_;
        dimensionedScalar tInf_;
        dimensionedScalar Cc_;
        dimensionedScalar Cv_;

        dimensionedScalar p0_;

        dimensionedScalar mcCoeff_;


public:

    //- Runtime type information
    TypeName("Merkle");


    // Constructors

        //- Construct for mixture
        Merkle
        (
            const dictionary& dict,
            const compressibleTwoPhaseMixture& mixture
        );


    //- Destructor
    virtual ~Merkle()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotP() const;

        //- Correct the Merkle phaseChange model
        virtual void correct();

        //- Read the dictionary and update
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cavitationModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
