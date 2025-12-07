import React from 'react';
import { LucideIcon } from 'lucide-react';

interface ResultCardProps {
  label: string;
  value: string | number;
  unit?: string;
  icon: LucideIcon;
  subValue?: string;
}

export const ResultCard: React.FC<ResultCardProps> = ({ label, value, unit, icon: Icon, subValue }) => {
  return (
    <div className="bg-slate-800/50 p-4 rounded-xl border border-slate-700 flex items-center space-x-4">
      <div className="p-3 bg-indigo-500/10 rounded-lg">
        <Icon className="w-6 h-6 text-indigo-400" />
      </div>
      <div>
        <p className="text-xs font-medium text-slate-400 uppercase">{label}</p>
        <p className="text-2xl font-bold text-white tracking-tight">
          {value} <span className="text-sm text-slate-500 font-normal">{unit}</span>
        </p>
        {subValue && <p className="text-xs text-slate-500 mt-1">{subValue}</p>}
      </div>
    </div>
  );
};