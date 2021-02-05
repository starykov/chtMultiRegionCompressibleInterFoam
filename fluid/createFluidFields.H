

PtrList<twoPhaseMixtureThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
// PtrList<volScalarField> rho1Fluid(fluidRegions.size());
// PtrList<volScalarField> rho2Fluid(fluidRegions.size());
// PtrList<volScalarField> alpha1Fluid(fluidRegions.size());
// PtrList<volScalarField> alpha2Fluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volScalarField> psi1Fluid(fluidRegions.size());
PtrList<volScalarField> psi2Fluid(fluidRegions.size());



PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());
PtrList<volScalarField> dgdtFluid(fluidRegions.size());


PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
//PtrList<uniformDimensionedScalarField> pRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressibleInterPhaseTransportModel>
    turbulenceFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());


PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

PtrList<bool> alphaRestartFluid(fluidRegions.size());
//PtrList<IOobject> alphaPhi10HeaderFluid(fluidRegions.size());
PtrList<surfaceScalarField> alphaPhi10Fluid(fluidRegions.size());


PtrList<tmp<surfaceScalarField>> talphaPhi1Corr0Fluid;

forAll(fluidRegions, i)
{
    Info<< "    Adding to pFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(UFluid[i])
        )
    );

    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;
    thermoFluid.set(i, new twoPhaseMixtureThermo(UFluid[i], phiFluid[i]));

    // Info<< "    Adding to alpha1Fluid\n" << endl;
    // alpha1Fluid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         "alpha1",
    //         thermoFluid[i].alpha1()
    //     )
    // );

    // Info<< "    Adding to alpha2Fluid\n" << endl;
    // alpha2Fluid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         "alpha1",
    //         thermoFluid[i].alpha2()
    //     )
    // );

    // Info<< "    Adding to rho1Fluid\n" << endl;
    // rho1Fluid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         "alpha1",
    //         thermoFluid[i].thermo1().rho()
    //     )
    // );

    // Info<< "    Adding to rho2Fluid\n" << endl;
    // rho2Fluid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         "alpha1",
    //         thermoFluid[i].thermo2().rho()
    //     )
    // );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].alpha1()*thermoFluid[i].thermo1().rho()+thermoFluid[i].alpha2()*thermoFluid[i].thermo2().rho()
        )
    );

    Info<< "    Adding to pMinFluid\n" << endl;
    pMinFluid.set
    (
        i,
        new dimensionedScalar
        (
            "pMin",
            dimPressure,
            thermoFluid[i]
        )
    );

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
    fluidRegions[i].setFluxRequired(thermoFluid[i].alpha1().name());

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero)
        )
    );

    Info<< "    Adding ghRef\n" << endl;
    dimensionedScalar ghRef(- mag(gFluid[i])*hRefFluid[i]);

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoFluid[i])*phiFluid[i]
        )
    );

    // Info<< "    Adding to dgdtFluid\n" << endl;
    // dgdtFluid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         "dgdt",
    //         thermoFluid[i].alpha1()*fvc::div(phiFluid[i])
    //     )
    // );


Info<< "    Adding to alphaPhi10Header\n" << endl;
    IOobject alphaPhi10Header
        (
        IOobject::groupName("alphaPhi0", thermoFluid[i].alpha1().group()),
        runTime.timeName(),
        fluidRegions[i],
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );


// Info<< "    Adding to alphaRestartFluid\n" << endl;
// //    const bool alphaRestart = alphaPhi10Header.typeHeaderOk<surfaceScalarField>();
//     alphaRestartFluid.set(i, alphaPhi10Header.typeHeaderOk<surfaceScalarField>());

// Info << i << "\n" << endl;
// Info << alphaRestartFluid[i] << "\n" << endl; 


// Info<< "    Adding to alphaRestartFluid2\n" << endl;
//     if (alphaRestartFluid[i])
//     {
//     Info << "Restarting alpha" << endl;
//     }

    Info<< "    Adding to alphaPhi10Fluid\n" << endl;
    alphaPhi10Fluid.set
    (
        i,
        new surfaceScalarField
        (
            alphaPhi10Header,
            phiFluid[i]*fvc::interpolate(thermoFluid[i].alpha1())
        )
    );

// MULES Correction
//tmp<surfaceScalarField> talphaPhi1Corr0;

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        new compressibleInterPhaseTransportModel
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            rhoPhiFluid[i],
            alphaPhi10Fluid[i],
            thermoFluid[i]
        )
    );

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

  

}

