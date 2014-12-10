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
    Foam::multiComponentMixture

Description
    Foam::multiComponentMixture

SourceFiles
    multiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentMixture_H
#define multiComponentMixture_H

#include "basicMultiComponentMixture.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multiComponentMixture
:
    public basicMultiComponentMixture
{
    // Private data

        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Construct as copy (not implemented)
        multiComponentMixture(const multiComponentMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database and mesh
        multiComponentMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const fvMesh&
        );

        //- Construct from dictionary and mesh
        multiComponentMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~multiComponentMixture()
    {}


    // Member functions

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
