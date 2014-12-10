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
    Foam::dragModels::segregated

Description
    Segregated drag model for use in regions with no obvious dispersed phase.

    Reference:
    \verbatim
        "Towards the Numerical Simulation of Multi-scale Two-phase Flows"
        H Marschall
        PhD Thesis, TU München, 2011
    \endverbatim

SourceFiles
    segregated.C

\*---------------------------------------------------------------------------*/

#ifndef segregated_H
#define segregated_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class segregated Declaration
\*---------------------------------------------------------------------------*/

class segregated
:
    public dragModel
{
private:

    // Private data

        //- M coefficient
        const dimensionedScalar m_;

        //- N coefficient
        const dimensionedScalar n_;


public:

    //- Runtime type information
    TypeName("segregated");


    // Constructors

        //- Construct from components
        segregated
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~segregated();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;

        //- The drag function used in the momentum equation
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
