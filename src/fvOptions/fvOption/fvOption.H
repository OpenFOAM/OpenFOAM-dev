/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::fv::option

Description
    Finite volume options abtract base class.  Provides a base set of controls,
    e.g.

        type            scalarExplicitSource    // source type
        active          on;                     // on/off switch

Note:
    On evaluation, source/sink options are to be added to the equation R.H.S.

SourceFiles
    fvOption.C
    fvOptionIO.C

\*---------------------------------------------------------------------------*/

#ifndef fvOption_H
#define fvOption_H

#include "fvMatricesFwd.H"
#include "volFieldsFwd.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class option Declaration
\*---------------------------------------------------------------------------*/

class option
{
protected:

    // Protected data

        //- Source name
        const word name_;

        //- Model type
        const word modelType_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Top level source dictionary
        dictionary dict_;

        //- Dictionary containing source coefficients
        dictionary coeffs_;

        //- Source active flag
        bool active_;

        //- Field names to apply source to - populated by derived models
        wordList fieldNames_;

        //- Applied flag list - corresponds to each fieldNames_ entry
        List<bool> applied_;


    // Protected functions

        //- Flag to bypass the apply flag list checking
        virtual bool alwaysApply() const;


public:

    //- Runtime type information
    TypeName("option");


     // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            option,
            dictionary,
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (name, modelType, dict, mesh)
        );


    // Constructors

        //- Construct from components
        option
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<option> clone() const
        {
            notImplemented("autoPtr<option> clone() const");
            return autoPtr<option>(NULL);
        }

        //- Return pointer to new fvOption object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh
            const fvMesh& mesh_;

            const word& name_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const word& name
            )
            :
                mesh_(mesh),
                name_(name)
            {}

            autoPtr<option> operator()(Istream& is) const
            {
                //const word name(is);
                const dictionary dict(is);

                return autoPtr<option>
                (
                    option::New(name_, dict, mesh_)
                );
            }
        };


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<option> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~option();


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return dictionary
            inline const dictionary& coeffs() const;

            //- Return const access to the source active flag
            inline bool active() const;

            //- Set the applied flag to true for field index fieldI
            inline void setApplied(const label fieldI);


        // Edit

            //- Return access to the source active flag
            inline bool& active();


        // Checks

            //- Is the source active?
            virtual bool isActive();

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;

            //- Check that the source has been applied
            virtual void checkApplied() const;


        // Evaluation

            // Correct

                //- Scalar
                virtual void correct(volScalarField& fld);

                //- Vector
                virtual void correct(volVectorField& fld);

                //- Spherical tensor
                virtual void correct(volSphericalTensorField& fld);

                //- Symmetric tensor
                virtual void correct(volSymmTensorField& fld);

                //- Tensor
                virtual void correct(volTensorField& fld);


            // Add explicit and implicit contributions

                //- Scalar
                virtual void addSup
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void addSup
                (
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void addSup
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void addSup
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void addSup
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


            // Add explicit and implicit contributions to compressible equations

                //- Scalar
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


            // Add explicit and implicit contributions to phase equations

                //- Scalar
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


            // Set values directly

                //- Scalar
                virtual void setValue
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void setValue
                (
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void setValue
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void setValue
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void setValue
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


        // IO

            //- Write the source header information
            virtual void writeHeader(Ostream&) const;

            //- Write the source footer information
            virtual void writeFooter(Ostream&) const;

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
