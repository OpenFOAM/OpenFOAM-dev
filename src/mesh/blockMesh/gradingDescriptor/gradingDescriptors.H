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
    Foam::gradingDescriptors

Description
    List of gradingDescriptor for the sections of a block with additional IO
    functionality

SourceFiles
    gradingDescriptors.C

\*---------------------------------------------------------------------------*/

#ifndef gradingDescriptors_H
#define gradingDescriptors_H

#include "gradingDescriptor.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class gradingDescriptors;
Istream& operator>>(Istream&, gradingDescriptors&);

/*---------------------------------------------------------------------------*\
                     Class gradingDescriptors Declaration
\*---------------------------------------------------------------------------*/

class gradingDescriptors
:
    public List<gradingDescriptor>
{

public:

    // Constructors

        //- Default constructor
        gradingDescriptors();

        //- Construct from a gradingDescriptor
        gradingDescriptors(const gradingDescriptor& gd);


    // Member functions

        //- Return the inverse gradingDescriptors with 1/expansionRatio
        gradingDescriptors inv() const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, gradingDescriptors&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
