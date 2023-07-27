/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::multicomponentMixture

Description
    Foam::multicomponentMixture

SourceFiles
    multicomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentMixture_H
#define multicomponentMixture_H

#include "PtrList.H"
#include "specieElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multicomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multicomponentMixture
{
public:

    // Public Typedefs

        //- The type of thermodynamics this mixture is instantiated for
        typedef ThermoType thermoType;


private:

    // Private Data

        //- List of specie thermo
        PtrList<ThermoType> specieThermos_;

        //- Table of species composition
        List<List<specieElement>> specieCompositions_;

        //- Cached IO error locations for delayed error messages
        List<IOerrorLocation> specieDictLocations_;


public:

    // Constructors

        //- Construct from a dictionary
        multicomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        multicomponentMixture
        (
            const multicomponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~multicomponentMixture()
    {}


    // Member Functions

        //- Return the specie names
        wordList specieNames() const;

        //- Read dictionary
        void read(const dictionary&);

        //- Return the raw specie thermodynamic data
        inline const PtrList<ThermoType>& specieThermos() const
        {
            return specieThermos_;
        }

        //- Return thermo based on index
        inline const ThermoType& specieThermo(const label speciei) const
        {
            return specieThermos_[speciei];
        }

        //- Return composition based on index
        const List<specieElement>& specieComposition(const label speciei) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multicomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
