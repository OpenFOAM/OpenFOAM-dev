/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::specie

Description
    Base class of the thermophysical property types.

SourceFiles
    specieI.H
    specie.C

\*---------------------------------------------------------------------------*/

#ifndef specie_H
#define specie_H

#include "word.H"
#include "scalar.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class specie Declaration
\*---------------------------------------------------------------------------*/

class specie
{
    // Private data

        //- Name of specie
        word name_;

        //- Number of moles of this component in the mixture
        scalar nMoles_;

        //- Molecular weight of specie [kg/kmol]
        scalar molWeight_;


public:

    //- Runtime type information
    ClassName("specie");


    // Public constants

        // Thermodynamic constants

            //- Universal gas constant [J/(kmol K)]
            static const scalar RR;

            //- Standard pressure [Pa]
            static const scalar Pstd;

            //- Standard temperature [K]
            static const scalar Tstd;


    // Constructors


        //- Construct from components without name
        inline specie(const scalar nMoles, const scalar molWeight);

        //- Construct from components with name
        inline specie
        (
            const word& name,
            const scalar nMoles,
            const scalar molWeight
        );

        //- Construct as copy
        inline specie(const specie&);

        //- Construct as named copy
        inline specie(const word& name, const specie&);

        //- Construct from Istream
        specie(Istream&);

        //- Construct from dictionary
        specie(const dictionary& dict);


    // Member Functions

        // Access

            //- Name
            inline const word& name() const;

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- No of moles of this species in mixture
            inline scalar nMoles() const;

            //- Gas constant [J/(kg K)]
            inline scalar R() const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator=(const specie&);

        inline void operator+=(const specie&);
        inline void operator-=(const specie&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend specie operator+(const specie&, const specie&);
        inline friend specie operator-(const specie&, const specie&);

        inline friend specie operator*(const scalar, const specie&);

        inline friend specie operator==(const specie&, const specie&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const specie&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
