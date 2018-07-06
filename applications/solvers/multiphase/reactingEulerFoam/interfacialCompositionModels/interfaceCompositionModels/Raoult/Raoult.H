/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::interfaceCompositionModels::Raoult

Description
    Raoult's law of ideal mixing. A separate composition model is given for
    each species. The composition of a species is equal to the value given by
    the model scaled by the species fraction in the bulk of the other phase.

SourceFiles
    Raoult.C

\*---------------------------------------------------------------------------*/

#ifndef Raoult_H
#define Raoult_H

#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                         Class Raoult Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class Raoult
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Non-vapour species fraction
        volScalarField YNonVapour_;

        //- Non-vapour species fraction derivative w.r.t. temperature
        volScalarField YNonVapourPrime_;

        //- Species' individual composition models
        HashTable<autoPtr<interfaceCompositionModel>> speciesModels_;


public:

    //- Runtime type information
    TypeName("Raoult");

    // Constructors

        //- Construct from components
        Raoult
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Raoult();


    // Member Functions

        //- Update the composition
        virtual void update(const volScalarField& Tf);

        //- The interface species fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- The interface species fraction derivative w.r.t. temperature
        virtual tmp<volScalarField> YfPrime
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceCompositionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Raoult.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
