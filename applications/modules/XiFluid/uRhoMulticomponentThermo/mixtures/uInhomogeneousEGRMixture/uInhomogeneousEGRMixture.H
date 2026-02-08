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
    Foam::uInhomogeneousEGRMixture

Description
    Unburnt gas combustion mixture for inhomogeneous lean, stoichiometric
    or rich conditions with exhaust gas recirculation

    in which the fuel mass-fraction \c fu and exhaust gas recirculation
    mass-fraction \c egr define the composition.

SourceFiles
    uInhomogeneousEGRMixture.C

\*---------------------------------------------------------------------------*/

#ifndef uInhomogeneousEGRMixture_H
#define uInhomogeneousEGRMixture_H

#include "uMixture.H"
#include "speciesTable.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class uInhomogeneousEGRMixture Declaration
\*---------------------------------------------------------------------------*/

class uInhomogeneousEGRMixture
:
    public uMixture
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
    ClassName("uInhomogeneousEGR");

    //- Specie indices
    enum species { FU, EGR };


    // Constructors

        //- Construct from a dictionary
        uInhomogeneousEGRMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        uInhomogeneousEGRMixture(const uInhomogeneousEGRMixture&) = delete;


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

        //- Return the stoichiometric air/fuel ratio
        scalar stoicRatio() const
        {
            return stoicRatio_;
        }

        //- Return the fuel-oxidant equivalence ratio
        scalar Phi(const scalarFieldListSlice& Yu) const;

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uInhomogeneousEGRMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
