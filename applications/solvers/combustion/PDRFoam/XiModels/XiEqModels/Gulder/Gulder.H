/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::XiEqModels::Gulder

Description
    Simple Gulder model for XiEq based on Gulders correlation
    with a linear correction function to give a plausible profile for XiEq.

SourceFiles
    Gulder.C

\*---------------------------------------------------------------------------*/

#ifndef Gulder_H
#define Gulder_H

#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class Gulder Declaration
\*---------------------------------------------------------------------------*/

class Gulder
:
    public XiEqModel
{
    // Private data

        //- Model constant
        scalar XiEqCoef_;

        //- Minimum laminar burning velocity
        const dimensionedScalar SuMin_;

        //- Schelkin effect Model constant
        scalar uPrimeCoef_;

        //- Use sub-grid Schelkin effect
        bool subGridSchelkin_;


    // Private Member Functions

        //- Disallow copy construct
        Gulder(const Gulder&);

        //- Disallow default bitwise assignment
        void operator=(const Gulder&);


public:

    //- Runtime type information
    TypeName("Gulder");


    // Constructors

        //- Construct from components
        Gulder
        (
            const dictionary& XiEqProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~Gulder();


    // Member Functions

        //- Return the flame-wrinking XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
