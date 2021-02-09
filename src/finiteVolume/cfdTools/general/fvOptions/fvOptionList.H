/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fv::optionList

Description
    List of finite volume options

SourceFile
    optionList.C

\*---------------------------------------------------------------------------*/

#ifndef fvOptionList_H
#define fvOptionList_H

#include "fvOption.H"
#include "PtrListDictionary.H"
#include "GeometricField.H"
#include "geometricOneField.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

namespace fv
{
    class optionList;
}

Ostream& operator<<(Ostream& os, const fv::optionList& options);

namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class optionList Declaration
\*---------------------------------------------------------------------------*/

class optionList
:
    public PtrListDictionary<option>
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Time index to check that all defined sources have been applied
        mutable label checkTimeIndex_;

        //- Sets of the fields that have had sources added by the options
        mutable PtrList<wordHashSet> addedToFields_;

        //- Sets of the fields that have been constrained by the options
        mutable PtrList<wordHashSet> constrainedFields_;

        //- Sets of the fields that have been corrected by the options
        mutable PtrList<wordHashSet> correctedFields_;


    // Protected Member Functions

        //- Read options dictionary
        bool readOptions(const dictionary& dict);

        //- Check that all options have been applied
        void checkApplied() const;

        //- Return source for equation with specified name and dimensions
        template<class Type, class ... AlphaRhoFieldTypes>
        tmp<fvMatrix<Type>> source
        (
            GeometricField<Type, fvPatchField, volMesh>& field,
            const word& fieldName,
            const dimensionSet& ds,
            const AlphaRhoFieldTypes& ... alphaRhos
        ) const;


private:

    // Private Member Functions

        //- Return the "options" sub-dictionary if present otherwise return dict
        const dictionary& optionsDict(const dictionary& dict) const;


public:

    //- Runtime type information
    TypeName("optionList");


    // Constructors

        //- Construct null
        optionList(const fvMesh& mesh);

        //- Construct from mesh and dictionary
        optionList(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        optionList(const optionList&) = delete;


    //- Destructor
    virtual ~optionList()
    {}


    // Member Functions

        //- Reset the source list
        void reset(const dictionary& dict);


        // Checks

            //- Return true if an option adds a source term to the given
            //  field's transport equation
            virtual bool addsToField(const word& fieldName) const;

            //- Return true if an option constrains the given field
            virtual bool constrainsField(const word& fieldName) const;

            //- Return true if an option corrects the given field
            virtual bool correctsField(const word& fieldName) const;


        // Sources

            //- Return source for an equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for an equation with a specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            ) const;

            //- Return source for a compressible equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for a compressible equation with a specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            ) const;

            //- Return source for a phase equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for a phase equation with a specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            ) const;

            //- Return source for a phase equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const geometricOneField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for a phase equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const geometricOneField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for a phase equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const geometricOneField& alpha,
                const geometricOneField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for an equation with a second time derivative
            template<class Type>
            tmp<fvMatrix<Type>> d2dt2
            (
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for an equation with a second time derivative
            template<class Type>
            tmp<fvMatrix<Type>> d2dt2
            (
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            ) const;


        // Constraints

            //- Apply constraints to an equation
            template<class Type>
            void constrain(fvMatrix<Type>& eqn) const;


        // Correction

            //- Apply correction to a field
            template<class Type>
            void correct
            (
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;


        // Mesh changes

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);

            //- Update for mesh motion
            virtual bool movePoints();


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

            //- Write data to Ostream
            virtual bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const optionList& options
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const optionList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvOptionListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
