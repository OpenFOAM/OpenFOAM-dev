/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::bHomogeneousMixture

Description
    Burnt gas mixture for homogeneous lean, stoichiometric or rich conditions

SourceFiles
    bHomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef bHomogeneousMixture_H
#define bHomogeneousMixture_H

#include "bMixture.H"
#include "speciesTable.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class bHomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

class bHomogeneousMixture
:
    public bMixture
{

public:

    //- Runtime type information
    ClassName("bHomogeneous");


    // Constructors

        //- Construct from a dictionary
        bHomogeneousMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        bHomogeneousMixture(const bHomogeneousMixture&) = delete;


    //- Destructor
    virtual ~bHomogeneousMixture()
    {}


    // Member Functions

        //- Return the mixture type name
        virtual const word& mixtureType() const
        {
            return typeName;
        }

        //- Return the table of species
        const speciesTable& species() const
        {
            return speciesTable::null;
        }

        //- Access the specie active flags
        const boolList& speciesActive() const
        {
            NotImplemented;
            return boolList::null();
        }

        //- Access the specie active flags
        boolList& speciesActive()
        {
            NotImplemented;
            return const_cast<boolList&>(boolList::null());
        }

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bHomogeneousMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
