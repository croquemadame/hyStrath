/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
    Foam::rampInletFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an oscillating condition in terms of
    amplitude and timeDuration.

        \f[
            x_p = (1 + a sin(2 \pi f t))x_{ref} + x_o
        \f]

    where

    \vartable
        x_p     | patch values
        x_{ref} | patch reference values
        x_o     | patch offset values
        a       | amplitude
        f       | timeDuration [1/s]
        t       | time [s]
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        refValue     | reference value         | yes         |
        offset       | offset value            | no          | 0.0
        amplitude    | oscillation amplitude   | yes         |
        timeDuration    | oscillation timeDuration   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            rampInlet;
        refValue        uniform 5.0;
        offset          0.0;
        amplitude       constant 0.5;
        timeDuration       constant 10;
    }
    \endverbatim

Note
    The amplitude and timeDuration entries are Function1 types, able to describe
    time varying functions.  The example above gives the usage for supplying
    constant values.

SeeAlso
    Foam::Function1

SourceFiles
    rampInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef rampInletFvPatchField_H
#define rampInletFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class rampInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class rampInletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Reference value
        Field<Type> refValue_;

        //- Offset
        Type offset_;

        //- Amplitude
        autoPtr<Function1<scalar> > amplitude_;

        //- timeDuration
        autoPtr<Function1<scalar> > timeDuration_;

        //- Current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Return current scale
        scalar currentScale() const;


public:

    //- Runtime type information
    TypeName("rampInlet");


    // Constructors

        //- Construct from patch and internal field
        rampInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rampInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rampInletFvPatchField
        //  onto a new patch
        rampInletFvPatchField
        (
            const rampInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rampInletFvPatchField
        (
            const rampInletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new rampInletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rampInletFvPatchField
        (
            const rampInletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new rampInletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            const Field<Type>& refValue() const
            {
                return refValue_;
            }

            //- Return reference to the ref value to allow adjustment
            Field<Type>& refValue()
            {
                return refValue_;
            }

            //- Return amplitude
            scalar amplitude() const
            {
                const scalar t = this->db().time().timeOutputValue();
                return amplitude_->value(t);
            }

            // Returning value() as reference is causing rval/lval errors; Removed for now
            //scalar& amplitude()
            //{
            //    const scalar t = this->db().time().timeOutputValue();
            //    return amplitude_->value(t);
            //}

            //- Return timeDuration
            scalar timeDuration() const
            {
                const scalar t = this->db().time().timeOutputValue();
                return timeDuration_->value(t);
            }

            //scalar& timeDuration()
            //{
            //    const scalar t = this->db().time().timeOutputValue();
            //    return timeDuration_->value(t);
            //}


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
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

#ifdef NoRepository
#   include "rampInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
