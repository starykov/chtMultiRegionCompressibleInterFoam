    const fvMesh& mesh = fluidRegions[i];


    twoPhaseMixtureThermo& mixture = thermoFluid[i];
    

    volScalarField& alpha1 = mixture.alpha1();
    volScalarField& alpha2 = mixture.alpha2();



    const volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& rho2 = mixture.thermo2().rho();

    volScalarField& rho = rhoFluid[i];



    volVectorField& U = UFluid[i];



    surfaceScalarField& phi = phiFluid[i];
    surfaceScalarField& rhoPhi = rhoPhiFluid[i];



    compressibleInterPhaseTransportModel& turbulence = turbulenceFluid[i];



    volScalarField& K = KFluid[i];
    //volScalarField& dgdt = dgdtFluid[i];

    volScalarField dgdt = alpha1*fvc::div(phi);

    volScalarField& p_rgh = p_rghFluid[i];



    const uniformDimensionedVectorField& g = gFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    pimpleNoLoopControl& pimple = pimples.pimple(i);



    surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

    IOobject alphaPhi10Header
        (
        IOobject::groupName("alphaPhi0", alpha1.group()),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );


    //const bool& alphaRestart = alphaRestartFluid[i];
    const bool& alphaRestart = alphaPhi10Header.typeHeaderOk<surfaceScalarField>();



    bool LTS = fv::localEulerDdt::enabled(mesh);



    tmp<surfaceScalarField> talphaPhi1Corr0;


        volScalarField& p = mixture.p();

    volScalarField& T = mixture.T();

    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    autoPtr<surfaceVectorField> Uf;



    const dimensionedScalar& pMin = pMinFluid[i];

    // This solver does not support moving mesh but it uses the pressure
    // equation of one which does, so we need a dummy face-momentum field
