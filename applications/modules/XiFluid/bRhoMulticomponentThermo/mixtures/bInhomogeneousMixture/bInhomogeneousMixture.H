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
    Foam::bInhomogeneousMixture

Description
    Burnt gas combustion mixture for inhomogeneous lean, stoichiometric
    or rich conditions

    in which the unburnt fuel mass-fraction \c ft (mixture fraction) defines the
    composition.

SourceFiles
    bInhomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef bInhomogeneousMixture_H
#define bInhomogeneousMixture_H

#include "bMixture.H"
#include "speciesTable.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bInhomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

class bInhomogeneousMixture
:
    public bMixture
{
protected:

    // Protected Data

        //- Table of specie names
        speciesTable species_;

        //- Stoichiometric air/fuel ratio
        scalar stoicRatio_;

        //- List of specie active flags
        boolList active_;


public:

    //- Runtime type information
    ClassName("bInhomogeneous");

    //- Specie indices
    enum species { FT };


    // Constructors

        //- Construct from a dictionary
        bInhomogeneousMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        bInhomogeneousMixture(const bInhomogeneousMixture&) = delete;


    //- Destructor
    virtual ~bInhomogeneousMixture()
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
            return species_;
        }

        //- Access the specie active flags
        const boolList& speciesActive() const
        {
            return active_;
        }

        //- Access the specie active flags
        boolList& speciesActive()
        {
            return active_;
        }

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bInhomogeneousMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
