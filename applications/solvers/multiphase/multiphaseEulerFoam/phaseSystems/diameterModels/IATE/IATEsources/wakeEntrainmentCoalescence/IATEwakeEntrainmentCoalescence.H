/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::diameterModels::IATEsources::wakeEntrainmentCoalescence

Description
    Bubble coalescence due to wake entrainment IATE source as defined in paper:

    Reference:
    \verbatim
        Ishii, M., Kim, S., & Kelly, J. (2005).
        Development of interfacial area transport equation.
        Nuclear Engineering and Technology, 37(6), 525-536.
    \endverbatim

SourceFiles
    wakeEntrainmentCoalescence.C

\*---------------------------------------------------------------------------*/

#ifndef IATEwakeEntrainmentCoalescence_H
#define IATEwakeEntrainmentCoalescence_H

#include "IATEsource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{

/*---------------------------------------------------------------------------*\
                 Class wakeEntrainmentCoalescence Declaration
\*---------------------------------------------------------------------------*/

class wakeEntrainmentCoalescence
:
    public IATEsource
{
    // Private Data

        dimensionedScalar Cwe_;


public:

    //- Runtime type information
    TypeName("wakeEntrainmentCoalescence");

    // Constructors

        wakeEntrainmentCoalescence
        (
            const IATE& iate,
            const dictionary& dict
        );


    //- Destructor
    virtual ~wakeEntrainmentCoalescence()
    {}


    // Member Functions

        virtual tmp<fvScalarMatrix> R
        (
            const volScalarField& alphai,
            volScalarField& kappai
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace IATEsources
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
