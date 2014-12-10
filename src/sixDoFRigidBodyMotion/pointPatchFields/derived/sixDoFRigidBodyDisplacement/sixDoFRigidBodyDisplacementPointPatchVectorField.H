/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::sixDoFRigidBodyDisplacementPointPatchVectorField

Description
    Foam::sixDoFRigidBodyDisplacementPointPatchVectorField

SourceFiles
    sixDoFRigidBodyDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyDisplacementPointPatchVectorField_H
#define sixDoFRigidBodyDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class sixDoFRigidBodyDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        //- Six dof motion object
        sixDoFRigidBodyMotion motion_;

        //- Initial positions of points on the patch
        pointField initialPoints_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rhoName == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- State of gravity lookup:
        //  -1 = not determined yet, as the BC may be instantiated before g has
        //       been read into the db yet.  Determination deferred until first
        //       call to updateCoeffs.  A g keyword was not supplied to the
        //       dictionary.
        //  -2 = as for -1, but a gravity value was specified in the dictionary,
        //       specifying a value in the dictionary is considered a fatal
        //       error if g is available from the db.
        //   0 = Use this boundary condition's own value of gravity, as not
        //       available from the db.
        //   1 = Lookup gravity from db.
        label lookupGravity_;

        //- Gravity vector to store when not available from the db
        vector g_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyDisplacement");


    // Constructors

        //- Construct from patch and internal field
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const sixDoFRigidBodyDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new sixDoFRigidBodyDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const sixDoFRigidBodyDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new sixDoFRigidBodyDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
