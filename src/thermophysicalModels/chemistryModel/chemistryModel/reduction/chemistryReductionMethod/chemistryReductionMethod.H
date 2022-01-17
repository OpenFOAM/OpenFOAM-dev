/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::chemistryReductionMethod

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    chemistryReductionMethod.C
    chemistryReductionMethods.C
    chemistryReductionMethodNew.C
    chemistryReductionMethodI.H

\*---------------------------------------------------------------------------*/

#ifndef chemistryReductionMethod_H
#define chemistryReductionMethod_H

#include "IOdictionary.H"
#include "DynamicField.H"
#include "Switch.H"
#include "cpuTime.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ThermoType>
class chemistryModel;

/*---------------------------------------------------------------------------*\
                   Class chemistryReductionMethod Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryReductionMethod
{
protected:

    // Protected Data

        //- Dictionary that store the algorithm data
        const dictionary coeffsDict_;

        //- Reference to the chemistry model
        chemistryModel<ThermoType>& chemistry_;

        //- Total number of species
        const label nSpecie_;

        //- Number of active species
        label nActiveSpecies_;

        //- List of disabled reactions (disabled = true)
        Field<bool> reactionsDisabled_;

        //- List of active species (active = true)
        List<bool> activeSpecies_;


    //- Protected Member Functions

        //- Initialise reduction of the mechanism
        void initReduceMechanism();

        //- End reduction of the mechanism
        void endReduceMechanism(List<label>& ctos, DynamicList<label>& stoc);


private:

    // Private Data

        //- Switch to select performance logging
        Switch log_;

        //- Tolerance for the mechanism reduction algorithm
        scalar tolerance_;

        //- CPU time for logging
        cpuTime cpuTime_;

        //- ...
        int64_t sumnActiveSpecies_;

        //- ...
        int64_t sumn_;

        //- ...
        scalar reduceMechCpuTime_;

        // Log file for the average time spent reducing the chemistry
        autoPtr<OFstream> cpuReduceFile_;

        // Write average number of species
        autoPtr<OFstream> nActiveSpeciesFile_;


public:

    //- Runtime type information
    TypeName("chemistryReductionMethod");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryReductionMethod,
        dictionary,
        (
            const IOdictionary& dict,
            chemistryModel<ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        chemistryReductionMethod
        (
            chemistryModel<ThermoType>& chemistry
        );

        //- Construct from components
        chemistryReductionMethod
        (
            const IOdictionary& dict,
            chemistryModel<ThermoType>& chemistry
        );


    // Selector

        static autoPtr<chemistryReductionMethod<ThermoType>> New
        (
            const IOdictionary& dict,
            chemistryModel<ThermoType>& chemistry
        );


    //- Destructor
    virtual ~chemistryReductionMethod();


    // Member Functions

        //- Return mechanism reduction active
        virtual bool active() const
        {
            return true;
        }

        //- Return the number of species
        inline label nSpecie();

        //- Return the number of active species
        inline label nActiveSpecies() const;

        //- Return the tolerance
        inline scalar tolerance() const;

        //- Return the active species
        inline const List<bool>& activeSpecies() const;

        //- Return whether or not a reaction is disabled
        inline bool reactionDisabled(const label i) const;

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            List<label>& ctos,
            DynamicList<label>& stoc,
            const label li
        ) = 0;

        //- ...
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryReductionMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryReductionMethod.C"
    #include "chemistryReductionMethodNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
