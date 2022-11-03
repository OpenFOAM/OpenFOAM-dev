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
    Foam::phaseInterface

Description
    Class to represent an interface between phases. Derivations can further
    specify the configuration of that interface; e.g., representing dispersal,
    displacement or sidedness.

SourceFiles
    phaseInterface.C
    phaseInterfaceI.H

\*---------------------------------------------------------------------------*/

#ifndef phaseInterface_H
#define phaseInterface_H

#include "phaseModel.H"
#include "uniformDimensionedFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseInterfaceKey;

/*---------------------------------------------------------------------------*\
                         Class phaseInterface Declaration
\*---------------------------------------------------------------------------*/

class phaseInterface
{
private:

    // Private Static Data

        //- List of head separators. This is the first separator in an
        //  interface name, and it delimits the names of phases on either side
        //  of the interface. Other separators are modifiers (e.g.,
        //  displacedBy, inThe, ...).
        static wordList headSeparators_;

        //- Map from old-format separators to separators
        static HashTable<word> oldSeparatorToSeparator_;


public:

    // Public Static Functions

        //- Get a reference to phase1 after sorting the phases by index
        static const phaseModel& getPhase1
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        );

        //- Get a reference to phase2 after sorting the phases by index
        static const phaseModel& getPhase2
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        );

        //- Add a head separator to the list
        static bool addHeadSeparator(const word& separator);

        //- Add a old separator to separator to the table
        static bool addOldSeparatorToSeparator
        (
            const word& oldSeparator,
            const word& separator
        );

        //- Split an interface name and return all its parts
        static wordList nameToNameParts
        (
            const phaseSystem& fluid,
            const word& name
        );

        //- Split an interface name and return its separators
        static wordList nameToSeparators
        (
            const phaseSystem& fluid,
            const word& name
        );

        //- Convert a list of separators into a type name
        static word separatorsToTypeName
        (
            const wordList& separators
        );

        //- Convert an interface name into a type name. Essentially just
        //  replaces valid phase names with a "<phase>" placeholder.
        static word nameToTypeName
        (
            const phaseSystem& fluid,
            const word& name
        );

        //- Convert interface name parts to an interface name
        static word namePartsToName
        (
            const phaseSystem& fluid,
            const wordList& nameParts
        );

        //- Convert old-format interface name parts to an interface name. Used
        //  in phaseSystem to provide backwards compatible input.
        static word oldNamePartsToName
        (
            const phaseSystem& fluid,
            const wordList& oldNameParts
        );

        //- Return references to the phases associated with a given name, and a
        //  list of valid separators
        static Tuple2<const phaseModel&, const phaseModel&> identifyPhases
        (
            const phaseSystem& fluid,
            const word& name,
            const wordList& separators
        );


private:

    // Private Data

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;


public:

    //- Runtime type information
    TypeName("phaseInterface");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseInterface,
            word,
            (
                const phaseSystem& fluid,
                const word& name
            ),
            (fluid, name)
        );


    // Constructors

        //- Construct from phases
        phaseInterface
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        );

        //- Construct from phases
        phaseInterface
        (
            const Tuple2<const phaseModel&, const phaseModel&>& phases
        );

        //- Construct from fluid and name
        phaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );

        //- Construct from fluid and key
        phaseInterface
        (
            const phaseSystem& fluid,
            const phaseInterfaceKey& name
        );

        //- Clone function
        virtual autoPtr<phaseInterface> clone() const;


    //- Destructor
    virtual ~phaseInterface();


    // Selectors

        //- Select given fluid and name
        static autoPtr<phaseInterface> New
        (
            const phaseSystem& fluid,
            const word& name
        );

        //- Select by combining two interfaces
        static autoPtr<phaseInterface> New
        (
            const phaseInterface& interface1,
            const phaseInterface& interface2
        );

        //- Class used for construction of PtrLists of phaseInterfaces
        class iNew
        {
            const phaseSystem& fluid_;

        public:

            iNew(const phaseSystem& fluid)
            :
                fluid_(fluid)
            {}

            autoPtr<phaseInterface> operator()(Istream& is) const
            {
                return phaseInterface::New(fluid_, word(is));
            }
        };


    // Static Member Functions

        //- Return the separator that delimits this interface's name
        static word separator()
        {
            return word::null;
        }


    // Member Functions

        //- Name
        virtual word name() const;

        //- Cast to derived type for use in a model
        template<class ModelType, class Derived>
        const Derived& modelCast() const
        {
            if (!isA<Derived>(*this))
            {
                FatalErrorInFunction
                    << "Constructing " << ModelType::typeName
                    << " for interface " << name()
                    << " which is not of the required type "
                    << Derived::typeName << exit(FatalError);
            }

            return refCast<const Derived>(*this);
        }


        // Access

            //- Return phase 1
            inline const phaseModel& phase1() const;

            //- Return phase 2
            inline const phaseModel& phase2() const;

            //- Return true if this phaseInterface contains the given phase
            inline bool contains(const phaseModel& phase) const;

            //- Return the other phase relative to the given phase
            //  Generates a FatalError if this phaseInterface does not contain
            //  the given phase
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return the index of the given phase. Generates a FatalError if
            //  this phaseInterface does not contain the given phase
            inline label index(const phaseModel& phase) const;

            //- Return the phase system
            inline const phaseSystem& fluid() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return gravitational acceleration
            inline const uniformDimensionedVectorField& g() const;


        // Properties

            //- Average density
            tmp<volScalarField> rho() const;

            //- Relative velocity magnitude
            tmp<volScalarField> magUr() const;

            //- Surface tension coefficient
            tmp<volScalarField> sigma() const;


        //- STL const_iterator
        class const_iterator
        {
            // Private Data

                //- Reference to the pair for which this is an iterator
                const phaseInterface& pair_;

                //- Current index
                label index_;

                //- Construct an iterator with the given index
                inline const_iterator(const phaseInterface&, const label index);

        public:

            friend class phaseInterface;

            // Constructors

                //- Construct from pair, moving to its 'begin' position
                inline explicit const_iterator(const phaseInterface&);


            // Access

                //- Return the current index
                inline label index() const;


            // Member Operators

                inline bool operator==(const const_iterator&) const;

                inline bool operator!=(const const_iterator&) const;

                inline const phaseModel& operator*() const;
                inline const phaseModel& operator()() const;

                inline const phaseModel& otherPhase() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);
        };


        //- const_iterator set to the beginning of the pair
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the pair
        inline const_iterator cend() const;

        //- const_iterator set to the beginning of the pair
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the pair
        inline const_iterator end() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseInterfaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
