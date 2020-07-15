/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::diameterModels::SecondaryPropertyModel

Description
    Base class for modeling the transport of an additional secondary property in
    between size groups by means of coalescence, breakup, drift and nucleation.
    The default behavior is conserved transport of the property, which can be
    overwritten in a derived class.

SourceFiles
    SecondaryPropertyModel.C

\*---------------------------------------------------------------------------*/

#ifndef SecondaryPropertyModel_H
#define SecondaryPropertyModel_H

#include "sizeGroup.H"
#include "driftModel.H"
#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                   Class SecondaryPropertyModel Declaration
\*---------------------------------------------------------------------------*/

template<class ModelType>
class SecondaryPropertyModel
:
    public ModelType,
    public regIOobject
{

protected:

    // Protected typedefs

        using SpTable = HashTable<const SecondaryPropertyModel<ModelType>*>;

    // Protected Data

        //- Reference to sizeGroup
        const sizeGroup& sizeGroup_;

        //- Table with pointers to all secondary properties of ModelType
        //  available in the registry
        SpTable SecondaryPropertyModelTable_;


public:

    // Constructors

        //- Construct from dictionary and sizeGroup
        SecondaryPropertyModel
        (
            const dictionary& dict,
            const sizeGroup& group
        );

        //- Disallow default bitwise copy construction
        SecondaryPropertyModel
        (
            const SecondaryPropertyModel<ModelType>&
        ) = delete;


    //- Destructor
    virtual ~SecondaryPropertyModel();


    // Member Functions

        // Access

            //- Return table with pointers to all secondary properties of
            //  ModelType available in the registry
            const SpTable& SecondaryPropertyModelTable();

            //- Return name of the secondary property of another class for
            //  for lookup in registry
            const word SecondaryPropertyName(const sizeGroup& fi) const;

            //- Return reference to secondary property field
            virtual const volScalarField& fld() const = 0;

            //- Access to secondary property source
            virtual volScalarField& src() = 0;

        // Edit

            //- Add coalescence contribution to secondary property source
            virtual void addCoalescence
            (
                const volScalarField& Su,
                const sizeGroup& fj,
                const sizeGroup& fk
            );

            //- Add breakup contribution to secondary property source
            virtual void addBreakup
            (
                const volScalarField& Su,
                const sizeGroup& fj
            );

            //- Add drift contribution to secondary property source
            virtual void addDrift
            (
                const volScalarField& Su,
                const sizeGroup& fu,
                const driftModel& model
            );

            //- Add drift contribution to secondary property source
            virtual void addNucleation
            (
                const volScalarField& Su,
                const sizeGroup& fi,
                const nucleationModel& model
            );

            //- Reset secondary property source
            virtual void reset();

        // Write

            //- Dummy write for regIOobject
            bool writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SecondaryPropertyModel<ModelType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SecondaryPropertyModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
