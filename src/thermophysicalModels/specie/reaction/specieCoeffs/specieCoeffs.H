/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
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
    Foam::specieCoeffs

Description
    Class to hold the specie index and its coefficients in the
    reaction rate expression.

SourceFiles
    specieCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef specieCoeffs_H
#define specieCoeffs_H

#include "specieExponent.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class specieCoeffs Declaration
\*---------------------------------------------------------------------------*/

class specieCoeffs
{
public:

    // Public Data

        //- Index of the specie
        label index;

        //- Stoichiometric coefficient
        scalar stoichCoeff;

        //- Exponent of the specie concentration
        specieExponent exponent;


    // Constructors

        //- Construct null
        specieCoeffs();

        //- Construct from species table and input stream
        specieCoeffs(const speciesTable& species, Istream& is);


    // Member Operators

        //- Equality comparison
        bool operator==(const specieCoeffs& sc) const
        {
            return index == sc.index;
        }

        //- Inequality comparison
        bool operator!=(const specieCoeffs& sc) const
        {
            return index != sc.index;
        }

        //- Write to output stream
        void write(const speciesTable& species, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
