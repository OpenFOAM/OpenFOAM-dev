/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::chemistryReduction

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    chemistryReduction.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryReduction_H
#define chemistryReduction_H

#include "chemistryReductionMethod.H"
#include "Switch.H"
#include "int64.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryReduction
:
    public chemistryReductionMethod<ThermoType>
{

protected:

    // Protected Data

        //- Dictionary that store the algorithm data
        const dictionary coeffsDict_;

        chemistryModel<ThermoType>& chemistry_;

        //- List of active species (active = true)
        List<bool> activeSpecies_;

        void initReduceMechanism();

        void endReduceMechanism();


private:

        //- Switch to select performance logging
        Switch log_;

        //- Tolerance for the mechanism reduction algorithm
        scalar tolerance_;

        const clockTime clockTime_;

        int64_t sumnActiveSpecies_;
        int64_t sumn_;
        scalar reduceMechCpuTime_;

        // Log file for the average time spent reducing the chemistry
        autoPtr<OFstream> cpuReduceFile_;

        // Write average number of species
        autoPtr<OFstream> nActiveSpeciesFile_;


public:

    // Constructors

        //- Construct from components
        chemistryReduction
        (
            const IOdictionary& dict,
            chemistryModel<ThermoType>& chemistry
        );


    //- Destructor
    virtual ~chemistryReduction();


    // Member Functions

        //- Return mechanism reduction active
        virtual bool active() const
        {
            return true;
        }

        //- Return the tolerance
        inline scalar tolerance() const;

        //- Return the active species
        inline const List<bool>& activeSpecies() const;

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            DynamicField<scalar>& sc,
            List<label>& ctos,
            DynamicList<label>& stoc,
            const label li
        ) = 0;

        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryReductionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryReduction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
