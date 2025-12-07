export enum ExposureOption {
  PORT_ONLY = 'Port Area Only',
  PORT_AND_ENDS = 'Port + Two Ends',
}

export interface MotorInputs {
  // Dimensions
  totalLength: number;        // mm
  chamberDiameter: number;    // mm
  propellantOD: number;       // mm
  propellantCore: number;     // mm (ID)
  propellantLength: number;   // mm
  exposureOption: ExposureOption;
  throatDiameter: number;     // mm
  throatErosion: number;      // mm (Total increase in diameter)
  
  // Nozzle Geometry
  nozzleConvergingAngle: number; // Beta (degrees)
  nozzleDivergingAngle: number;  // Alpha (degrees)
  nozzleExpansionRatio: number;  // Ae / At
  nozzleEfficiency: number;      // dimensionless (0-1 typically)

  // Properties
  density: number;            // g/cm^3 (gm/cc)
  specificHeatRatio: number;  // k (dimensionless)
  temperature: number;        // K
  molecularMass: number;      // g/mol
  propellantGStar: number;    // Propellant Threshold
  atmosphericPressure: number; // MPa
  
  // Burn Rate Law: r = a * P^n
  burnRateA: number;
  burnRateN: number;
}

export interface SimulationStep {
  step: number;
  webRegression: number;      // mm
  webRemaining: number;       // mm
  kn: number;                 // dimensionless
  burningArea: number;        // mm^2
  throatDiameter: number;     // mm
  pressure: number;           // MPa
  thrust: number;             // N
  time: number;               // s
}

export interface BurnRateStep {
  pressure: number;           // MPa
  rate: number;               // mm/s
}

export interface SimulationResults {
  chamberVolume: number;      // mm^3
  propellantVolume: number;   // mm^3
  propellantMass: number;     // g
  initialKn: number;
  maxKn: number;
  maxPressure: number;        // MPa
  maxThrust: number;          // N
  averageThrust: number;      // N (during burn)
  specificImpulse: number;    // s
  burnTime: number;           // s
  burnOutTime: number;        // s
  webThickness: number;       // mm
  initialBurningArea: number; // mm^2
  finalThroatDiameter: number;// mm
  nozzleExitDiameter: number; // mm
  gasConstantR: number;       // J/(kg*K)
  cStar: number;              // m/s
  data: SimulationStep[];
  burnRateData: BurnRateStep[];
}