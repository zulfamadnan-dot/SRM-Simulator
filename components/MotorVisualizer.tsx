import React, { useMemo } from 'react';
import { MotorInputs } from '../types';

interface MotorVisualizerProps {
  inputs: MotorInputs;
  calculatedExitDiameter: number;
}

export const MotorVisualizer: React.FC<MotorVisualizerProps> = ({ inputs, calculatedExitDiameter }) => {
  const {
    totalLength, // Chamber length
    chamberDiameter,
    propellantLength,
    propellantOD,
    propellantCore,
    throatDiameter,
    nozzleConvergingAngle,
    nozzleDivergingAngle,
  } = inputs;

  // Geometry Calculations
  const geometry = useMemo(() => {
    // Prevent division by zero or extreme values
    const beta = Math.max(5, Math.min(89, nozzleConvergingAngle)) * (Math.PI / 180);
    const alpha = Math.max(1, Math.min(89, nozzleDivergingAngle)) * (Math.PI / 180);

    const rChamber = chamberDiameter / 2;
    const rThroat = throatDiameter / 2;
    const rExit = calculatedExitDiameter / 2;
    const rPropOuter = propellantOD / 2;
    const rPropCore = propellantCore / 2;

    // Lengths
    const lConv = Math.abs((rChamber - rThroat) / Math.tan(beta));
    const lDiv = Math.abs((rExit - rThroat) / Math.tan(alpha));
    const lNozzle = lConv + lDiv;
    const lTotalVisual = totalLength + lNozzle;
    
    // Centering Offset for Propellant
    const grainOffset = (totalLength - propellantLength) / 2;

    // Viewbox sizing
    const maxDiameter = Math.max(chamberDiameter, calculatedExitDiameter);
    const paddingX = 40; // mm
    const paddingY = 60; // mm
    const viewBoxHeight = maxDiameter + paddingY * 2;
    const viewBoxWidth = lTotalVisual + paddingX * 2;
    
    const centerY = viewBoxHeight / 2;
    const startX = paddingX;

    return {
        rChamber, rThroat, rExit, rPropOuter, rPropCore,
        lConv, lDiv, lNozzle, lTotalVisual, grainOffset,
        viewBoxWidth, viewBoxHeight, centerY, startX
    };
  }, [inputs, calculatedExitDiameter]);

  const {
      rChamber, rThroat, rExit, rPropOuter, rPropCore,
      lConv, lDiv, lTotalVisual, grainOffset,
      viewBoxWidth, viewBoxHeight, centerY, startX
  } = geometry;

  // Helper for formatting numbers
  const fmt = (n: number) => n.toFixed(1);

  return (
    <div className="w-full bg-slate-900 border border-slate-800 rounded-xl p-4 shadow-xl overflow-hidden flex flex-col items-center">
      <h3 className="text-slate-400 text-sm font-semibold uppercase tracking-wider mb-4 self-start">Motor Configuration</h3>
      <div className="w-full overflow-x-auto">
      <svg 
        width="100%" 
        height="300px" 
        viewBox={`0 0 ${viewBoxWidth} ${viewBoxHeight}`} 
        preserveAspectRatio="xMidYMid meet"
        className="mx-auto"
      >
        <defs>
            <pattern id="propellantPattern" patternUnits="userSpaceOnUse" width="4" height="4" patternTransform="rotate(45)">
                <rect width="2" height="4" transform="translate(0,0)" fill="#10b981" opacity="0.2" />
            </pattern>
            <marker id="arrow-start" markerWidth="10" markerHeight="10" refX="0" refY="3" orient="auto" markerUnits="strokeWidth">
              <path d="M9,0 L0,3 L9,6 L9,0" fill="#64748b" />
            </marker>
            <marker id="arrow-end" markerWidth="10" markerHeight="10" refX="10" refY="3" orient="auto" markerUnits="strokeWidth">
              <path d="M0,0 L10,3 L0,6 L0,0" fill="#64748b" />
            </marker>
        </defs>

        {/* Centerline */}
        <line 
            x1={0} 
            y1={centerY} 
            x2={viewBoxWidth} 
            y2={centerY} 
            stroke="#334155" 
            strokeDasharray="10,5" 
            strokeWidth="1" 
        />

        {/* --- PROPELLANT --- */}
        <g transform={`translate(${startX + grainOffset}, ${centerY})`}>
            {/* Top Grain */}
            <rect 
                x={0} 
                y={-rPropOuter} 
                width={propellantLength} 
                height={rPropOuter - rPropCore} 
                fill="url(#propellantPattern)" 
                stroke="#059669" 
                strokeWidth="1"
            />
            {/* Bottom Grain */}
             <rect 
                x={0} 
                y={rPropCore} 
                width={propellantLength} 
                height={rPropOuter - rPropCore} 
                fill="url(#propellantPattern)" 
                stroke="#059669" 
                strokeWidth="1"
            />
        </g>

        {/* --- CHAMBER CASING --- */}
        <g transform={`translate(${startX}, ${centerY})`}>
             {/* Top Wall */}
             <line x1={0} y1={-rChamber} x2={totalLength} y2={-rChamber} stroke="#94a3b8" strokeWidth="2" />
             {/* Bottom Wall */}
             <line x1={0} y1={rChamber} x2={totalLength} y2={rChamber} stroke="#94a3b8" strokeWidth="2" />
             {/* Head End Closure */}
             <line x1={0} y1={-rChamber} x2={0} y2={rChamber} stroke="#94a3b8" strokeWidth="2" />
        </g>

        {/* --- NOZZLE --- */}
        <g transform={`translate(${startX + totalLength}, ${centerY})`}>
            {/* Top Path */}
            <path 
                d={`
                    M 0 ${-rChamber} 
                    L ${lConv} ${-rThroat} 
                    L ${lConv + lDiv} ${-rExit}
                `} 
                stroke="#f59e0b" 
                strokeWidth="2" 
                fill="none" 
            />
            {/* Bottom Path */}
             <path 
                d={`
                    M 0 ${rChamber} 
                    L ${lConv} ${rThroat} 
                    L ${lConv + lDiv} ${rExit}
                `} 
                stroke="#f59e0b" 
                strokeWidth="2" 
                fill="none" 
            />
            
            {/* Throat Line (dashed visual aid) */}
            <line x1={lConv} y1={-rThroat} x2={lConv} y2={rThroat} stroke="#f59e0b" strokeWidth="1" strokeDasharray="2,2" opacity="0.5" />
        </g>

        {/* --- DIMENSIONS & ANNOTATIONS --- */}
        
        {/* Total Chamber Length Dimension */}
        <g transform={`translate(${startX}, ${centerY + rChamber + 20})`}>
            <line x1={0} y1={0} x2={totalLength} y2={0} stroke="#64748b" strokeWidth="1" markerStart="url(#arrow-start)" markerEnd="url(#arrow-end)" />
            <text x={totalLength / 2} y={15} fill="#94a3b8" fontSize="10" textAnchor="middle">L_motor: {fmt(totalLength)}</text>
            {/* Drop lines */}
            <line x1={0} y1={-15} x2={0} y2={0} stroke="#64748b" strokeWidth="0.5" opacity="0.5" />
            <line x1={totalLength} y1={-15} x2={totalLength} y2={0} stroke="#64748b" strokeWidth="0.5" opacity="0.5" />
        </g>

        {/* Propellant Length Dimension */}
        <g transform={`translate(${startX + grainOffset}, ${centerY - rChamber - 20})`}>
            <line x1={0} y1={0} x2={propellantLength} y2={0} stroke="#64748b" strokeWidth="1" markerStart="url(#arrow-start)" markerEnd="url(#arrow-end)" />
            <text x={propellantLength / 2} y={-5} fill="#10b981" fontSize="10" textAnchor="middle">L_grain: {fmt(propellantLength)}</text>
             {/* Drop lines */}
            <line x1={0} y1={0} x2={0} y2={15} stroke="#64748b" strokeWidth="0.5" opacity="0.5" />
            <line x1={propellantLength} y1={0} x2={propellantLength} y2={15} stroke="#64748b" strokeWidth="0.5" opacity="0.5" />
        </g>
        
        {/* Throat Label */}
        <g transform={`translate(${startX + totalLength + lConv}, ${centerY})`}>
             <text x={0} y={rExit + 40} fill="#f59e0b" fontSize="10" textAnchor="middle">Throat: {fmt(throatDiameter)}</text>
             <line x1={0} y1={rThroat} x2={0} y2={rExit + 30} stroke="#f59e0b" strokeWidth="0.5" strokeDasharray="2,2" />
        </g>

        {/* Exit Label */}
        <g transform={`translate(${startX + totalLength + lConv + lDiv}, ${centerY})`}>
             <text x={0} y={rExit + 20} fill="#f59e0b" fontSize="10" textAnchor="middle">Exit: {fmt(calculatedExitDiameter)}</text>
             <line x1={0} y1={rExit} x2={0} y2={rExit + 10} stroke="#f59e0b" strokeWidth="0.5" />
        </g>

        {/* Nozzle Lengths */}
        <g transform={`translate(${startX + totalLength}, ${centerY - Math.max(rChamber, rExit) - 10})`}>
             <text x={lConv/2} y={0} fill="#64748b" fontSize="9" textAnchor="middle">Conv</text>
             <text x={lConv + lDiv/2} y={0} fill="#64748b" fontSize="9" textAnchor="middle">Div</text>
        </g>

      </svg>
      </div>
    </div>
  );
};