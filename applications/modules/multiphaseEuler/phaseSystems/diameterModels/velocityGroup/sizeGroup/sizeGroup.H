/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Foam::diameterModels::sizeGroup

Description
    Single size class fraction field representing a fixed particle volume as
    defined by the user through the corresponding sphere equivalent diameter.

    If present, the field is read from the start time directory, e.g. during
    restart. Otherwise, it is constructed as a uniform field using the size
    class fraction value provided by the user. In the latter case, the boundary
    condition types are taken from the 'f.<phaseName>' field, which represents
    the sum of all size group fractions of a phase. The user specified value is
    also applied at fixed value boundary conditions and used as inlet value for
    mixed boundary conditions.

    An alternative diameter field is provided by the selected shape model, e.g.
    a collisional diameter, which is then utilised in selected population
    balance submodels, e.g. for modelling fractal aggregation.

Usage
    \table
        Property     | Description                | Required    | Default value
        dSph         | Sphere equivalent diameter | yes         | none
        value        | Initial and boundary condition value of\\
                       size class fraction        | yes         | none
        shapeModel   | Shape model providing an alternative diameter field\\
                                                  | yes         | none
    \endtable

See also
    Foam::diameterModels::populationBalanceModel
    Foam::diameterModels::velocityGroup
    Foam::diameterModels::shapeModel

SourceFiles
    sizeGroup.C

\*---------------------------------------------------------------------------*/

#ifndef sizeGroup_H
#define sizeGroup_H

#include "dictionaryEntry.H"
#include "velocityGroup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

class shapeModel;

/*---------------------------------------------------------------------------*\
                          Class sizeGroup Declaration
\*---------------------------------------------------------------------------*/

class sizeGroup
:
    public volScalarField
{
private:

    // Private Data

        //- Label of this sizeGroup within the corresponding populationBalance
        const label i_;

        //- VelocityGroup this sizeGroup belongs to
        const velocityGroup& group_;

        //- Sphere equivalent diameter of the sizeGroup
        const dimensionedScalar dSph_;

        //- Representative volume of the sizeGroup
        const dimensionedScalar x_;

        //- Model for describing the representative shape of the elements in the
        //  sizeGroup
        autoPtr<shapeModel> shapeModel_;


public:

    // Static Member Functions

        //- Return IO for a size-group field. Construction helper.
        static IOobject fieldIo
        (
            const word& name,
            const label i,
            const velocityGroup& group,
            const IOobject::readOption r = IOobject::NO_READ,
            const bool registerObject = true
        );

        //- Return IO a size-group field. Construction helper.
        static tmp<volScalarField> field
        (
            const word& name,
            const label i,
            const velocityGroup& group
        );


    // Constructors

        //- Construct from index, dictionary and velocity group
        sizeGroup
        (
            const label i,
            const dictionary& dict,
            const velocityGroup& group
        );

        //- Return clone
        autoPtr<sizeGroup> clone() const;

        //- Return a pointer to a new sizeGroup created from Istream
        class iNew
        {
            const velocityGroup& group_;

            mutable label i_;


        public:

            iNew(const velocityGroup& group, const label i0)
            :
                group_(group),
                i_(i0)
            {}

            autoPtr<sizeGroup> operator()(Istream& is) const
            {
                token t(is);

                // Check and filter for old syntax (remove in due course)
                if (t.isWord())
                {
                    const word fName = "f" + Foam::name(i_);

                    if (t.wordToken() != fName)
                    {
                        FatalErrorInFunction
                            << "The name '" << t.wordToken() << "' should not "
                            << "have been provided for the fraction field "
                            << IOobject::groupName(fName, group_.phase().name())
                            << " of size-group #" << i_
                            << " of population balance " << group_.popBalName()
                            << ". Only the coefficients are required."
                            << exit(FatalError);
                    }
                }
                else
                {
                    is.putBack(t);
                }

                autoPtr<sizeGroup> result(new sizeGroup(i_, is, group_));

                i_ ++;

                return result;
            }
        };


    //- Destructor
    virtual ~sizeGroup();


    // Member Functions

        //- Return index of the size group within the population balance
        label i() const;

        //- Return const-reference to the phase
        inline const phaseModel& phase() const;

        //- Return const-reference to the velocityGroup
        inline const velocityGroup& group() const;

        //- Return representative spherical diameter of the sizeGroup
        inline const dimensionedScalar& dSph() const;

        //- Return representative volume of the sizeGroup
        inline const dimensionedScalar& x() const;

        //- Return reference to diameterModel of the phase
        inline autoPtr<shapeModel>& shapeModelPtr();

        //- Return const-reference to diameterModel of the phase
        inline const autoPtr<shapeModel>& shapeModelPtr() const;

        //- Return representative surface area of the sizeGroup
        const tmp<volScalarField> a() const;

        //- Return representative diameter of the sizeGroup
        const tmp<volScalarField> d() const;

        //- Correct secondary properties
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sizeGroupI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
