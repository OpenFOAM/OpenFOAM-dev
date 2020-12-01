/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Finite volume options abstract base class.  Provides a base set of
    controls, e.g.:
    \verbatim
        type            scalarExplicitSource    // source type
    \endverbatim

Note
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
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class mapPolyMesh;

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

        //- Field names to apply source to - populated by derived models
        wordList fieldNames_;

        //- Applied flag list - corresponds to each fieldNames_ entry
        mutable List<bool> applied_;


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
            NotImplemented;
            return autoPtr<option>(nullptr);
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
                // const word name(is);
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

            //- Set the applied flag to true for field index fieldi
            inline void setApplied(const label fieldi) const;


        // Checks

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;

            //- Check that the source has been applied
            virtual void checkApplied() const;


        // Evaluation

            // Explicit and implicit sources

                virtual void addSup
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                ) const;


            // Explicit and implicit sources for compressible equations

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                ) const;


            // Explicit and implicit sources for phase equations

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                ) const;


            // Constraints

                virtual void constrain
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                ) const;

                virtual void constrain
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                ) const;

                virtual void constrain
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void constrain
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                ) const;

                virtual void constrain
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                ) const;


            // Correction

                virtual void correct(volScalarField& field) const;
                virtual void correct(volVectorField& field) const;
                virtual void correct(volSphericalTensorField& field) const;
                virtual void correct(volSymmTensorField& field) const;
                virtual void correct(volTensorField& field) const;


        // Mesh changes

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);

            //- Update for mesh motion
            virtual bool movePoints();


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
