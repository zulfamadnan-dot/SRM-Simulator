import React from 'react';

interface InputGroupProps {
  label: string;
  name: string;
  value: number | string;
  onChange: (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => void;
  type?: 'number' | 'select';
  options?: string[];
  step?: string;
  unit?: string;
  description?: string;
}

export const InputGroup: React.FC<InputGroupProps> = ({
  label,
  name,
  value,
  onChange,
  type = 'number',
  options = [],
  step = "0.01",
  unit,
  description
}) => {
  return (
    <div className="flex flex-col space-y-1">
      <label className="text-xs font-medium text-slate-400 uppercase tracking-wider flex justify-between">
        {label}
        {unit && <span className="text-slate-500 normal-case lowercase">({unit})</span>}
      </label>
      {type === 'select' ? (
        <select
          name={name}
          value={value}
          onChange={onChange}
          className="bg-slate-900 border border-slate-700 text-slate-100 text-sm rounded-md focus:ring-indigo-500 focus:border-indigo-500 block w-full p-2.5 transition-colors"
        >
          {options.map((opt) => (
            <option key={opt} value={opt}>
              {opt}
            </option>
          ))}
        </select>
      ) : (
        <input
          type="number"
          name={name}
          value={value}
          onChange={onChange}
          step={step}
          className="bg-slate-900 border border-slate-700 text-slate-100 text-sm rounded-md focus:ring-indigo-500 focus:border-indigo-500 block w-full p-2.5 transition-colors"
        />
      )}
      {description && <p className="text-[10px] text-slate-600">{description}</p>}
    </div>
  );
};