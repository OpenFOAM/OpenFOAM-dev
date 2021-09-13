/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef chemistryReductionMethod_H
#define chemistryReductionMethod_H

#include "IOdictionary.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ThermoType>
class chemistryModel;

/*---------------------------------------------------------------------------*\
                       Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryReductionMethod
{
    // Private Data

        //- Cache the number of species
        const label nSpecie_;


protected:

    // Protected Data

        //- Number of active species
        label nActiveSpecies_;

        Field<bool> reactionsDisabled_;


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

        //- Is mechanism reduction active?
        virtual bool active() const = 0;

        //- Return the number of species
        inline label nSpecie();

        //- Return the number of active species
        inline label nActiveSpecies() const;

        //- Return the tolerance
        inline scalar tolerance() const;

        inline bool reactionDisabled(const label i) const
        {
            return reactionsDisabled_[i];
        }

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

        virtual void update() = 0;
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
