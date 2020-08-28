/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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
    Foam::dragModels::TomiyamaCorrelated

Description
    Correlation of Tomiyama et al.

    Reference:
    \verbatim
        Tomiyama, A., Celata, G. P., Hosokawa, S., & Yoshida, S. (2002).
        Terminal velocity of single bubbles in surface tension force dominant
        regime.
        International Journal of Multiphase Flow, 28(9), 1497-1519.
    \endverbatim

SourceFiles
    TomiyamaCorrelated.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaCorrelated_H
#define TomiyamaCorrelated_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                     Class TomiyamaCorrelated Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaCorrelated
:
    public dragModel
{
    // Private Data

        //- Coefficient
        const dimensionedScalar A_;


public:

    //- Runtime type information
    TypeName("TomiyamaCorrelated");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaCorrelated
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~TomiyamaCorrelated();


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
