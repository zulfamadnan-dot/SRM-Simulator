import { ExposureOption, MotorInputs } from './types';

export const INITIAL_INPUTS: MotorInputs = {
  totalLength: 300,
  chamberDiameter: 54,
  propellantOD: 50,
  propellantCore: 15,
  propellantLength: 250,
  exposureOption: ExposureOption.PORT_AND_ENDS,
  throatDiameter: 10,
  throatErosion: 1.0, // Default 1mm total erosion
  
  nozzleConvergingAngle: 30, // Beta
  nozzleDivergingAngle: 12,  // Alpha
  nozzleExpansionRatio: 6.0,
  nozzleEfficiency: 0.9,
  
  density: 1.7, // g/cm^3
  specificHeatRatio: 1.2,
  temperature: 3000,
  molecularMass: 24,
  propellantGStar: 150, // Default threshold
  atmosphericPressure: 0.101325, // MPa
  burnRateA: 0.5,
  burnRateN: 0.3,
};
