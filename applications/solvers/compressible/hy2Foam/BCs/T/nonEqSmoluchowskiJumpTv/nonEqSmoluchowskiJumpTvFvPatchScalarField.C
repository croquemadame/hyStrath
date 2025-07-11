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

\*---------------------------------------------------------------------------*/

#include "nonEqSmoluchowskiJumpTvFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mathematicalConstants.H"

#include <string.H>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::
nonEqSmoluchowskiJumpTvFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho"),
    muName_("mu"),
    //alphaName_("alphave"),
    gammaName_("gammatr"),
    mfpName_("mfp"),
    accommodationCoeff_(1.0),
    Twall_(p.size(), 0.0)
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
    alphaName_ = "alphave_" + specieName_;

    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::
nonEqSmoluchowskiJumpTvFvPatchScalarField
(
    const nonEqSmoluchowskiJumpTvFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    specieName_(ptf.specieName_),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_),
    muName_(ptf.muName_),
    alphaName_(ptf.alphaName_),
    gammaName_(ptf.gammaName_),
    mfpName_(ptf.mfpName_),
    accommodationCoeff_(ptf.accommodationCoeff_),
    Twall_(ptf.Twall_)
{}


Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::
nonEqSmoluchowskiJumpTvFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    muName_(dict.lookupOrDefault<word>("mu", "mu")),
    //alphaName_(dict.lookupOrDefault<word>("alphave", "alphave")),
    gammaName_(dict.lookupOrDefault<word>("gammatr", "gammatr")),
    mfpName_(dict.lookupOrDefault<word>("mfp", "mfp")),
    accommodationCoeff_(readScalar(dict.lookup("accommodationCoeff"))),
    Twall_("Twall", dict, p.size())
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
    alphaName_ = "alphave_" + specieName_;

    if
    (
        mag(accommodationCoeff_) < SMALL
     || mag(accommodationCoeff_) > 1.0
    )
    {
        FatalIOErrorIn
        (
            "nonEqSmoluchowskiJumpTvFvPatchScalarField::"
            "nonEqSmoluchowskiJumpTvFvPatchScalarField"
            "("
            "    const fvPatch&,"
            "    const DimensionedField<scalar, volMesh>&,"
            "    const dictionary&"
            ")",
            dict
        )   << "unphysical accommodationCoeff specified"
            << "(0 < accommodationCoeff <= 1)" << endl
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::
nonEqSmoluchowskiJumpTvFvPatchScalarField
(
    const nonEqSmoluchowskiJumpTvFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    accommodationCoeff_(ptpsf.accommodationCoeff_),
    Twall_(ptpsf.Twall_)
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& pmu =
        patch().lookupPatchField<volScalarField, scalar>(muName_);
    const fvPatchScalarField& palpha =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);
    const fvPatchScalarField& pgammatr =
        patch().lookupPatchField<volScalarField, scalar>(gammaName_);
    const fvPatchScalarField& pmfp =
        patch().lookupPatchField<volScalarField, scalar>(mfpName_);
    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    const fvPatchVectorField& pU =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    scalarField C2
    (
        pmfp*2.0*palpha/(pgammatr + 1.0)/(pmu/palpha)
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    );

//    scalarField aCoeff(prho.snGrad() - prho/C2);
//    scalarField KEbyRho(0.5*magSqr(pU));

    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C2));
    refValue() = Twall_;
    refGrad() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void Foam::nonEqSmoluchowskiJumpTvFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("mu", "mu", muName_);
    os.writeEntryIfDifferent<word>("alphave", "alphave", alphaName_);
    os.writeEntryIfDifferent<word>("gammatr", "gammatr", gammaName_);
    os.writeEntryIfDifferent<word>("mfp", "mfp", mfpName_);

    os.writeKeyword("accommodationCoeff")
        << accommodationCoeff_ << token::END_STATEMENT << nl;
    Twall_.writeEntry("Twall", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nonEqSmoluchowskiJumpTvFvPatchScalarField
    );
}


// ************************************************************************* //
