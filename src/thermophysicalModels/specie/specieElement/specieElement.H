/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::specieElement

Description

SourceFiles
    specieElementI.H

\*---------------------------------------------------------------------------*/

#ifndef specieElement_H
#define specieElement_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class specieElement;
Ostream& operator<<(Ostream&, const specieElement&);


/*---------------------------------------------------------------------------*\
                        Class specieElement Declaration
\*---------------------------------------------------------------------------*/

class specieElement
{
    // Private Data

        //- Name of the element
        word name_;

        //- Number of atoms of this element in the specie
        label nAtoms_;


public:

    // Constructors

        //- Construct null
        inline specieElement();

        //- Construct from components
        inline specieElement(const word& name, const label nAtoms);

        //- Construct from Istream
        inline specieElement(Istream&);


    // Member Functions

        //- Return the name of the element
        inline const word& name() const;

        //- Return non-const access to the name of the element
        inline word& name();

        //- Return the number of atoms of this element in the specie
        inline label nAtoms() const;

        //- Return non-const access to the number of atoms of this element
        //  in the specie
        inline label& nAtoms();


    // Member Operators

        //- Equality comparison
        inline bool operator==(const specieElement&) const;

        //- Inequality comparison
        inline bool operator!=(const specieElement&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const specieElement&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specieElementI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
