/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::diameterModels::IATEsources::randomCoalescence

Description
    Random coalescence IATE source as defined in paper:

    \verbatim
        "Development of Interfacial Area Transport Equation"
        M. Ishii,
        S. Kim,
        J Kelly,
        Nuclear Engineering and Technology, Vol.37 No.6 December 2005
    \endverbatim


SourceFiles
    randomCoalescence.C

\*---------------------------------------------------------------------------*/

#ifndef randomCoalescence_H
#define randomCoalescence_H

#include "IATEsource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{

/*---------------------------------------------------------------------------*\
                           Class randomCoalescence Declaration
\*---------------------------------------------------------------------------*/

class randomCoalescence
:
    public IATEsource
{
    // Private data

        dimensionedScalar Crc_;
        dimensionedScalar C_;
        dimensionedScalar alphaMax_;


public:

    //- Runtime type information
    TypeName("randomCoalescence");

    // Constructors

        randomCoalescence
        (
            const IATE& iate,
            const dictionary& dict
        );


    //- Destructor
    virtual ~randomCoalescence()
    {}


    // Member Functions

        virtual tmp<volScalarField> R() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace IATEsources
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
