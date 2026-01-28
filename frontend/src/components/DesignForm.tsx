import React, { useState } from 'react';

interface DesignFormProps {
  onSubmit: (data: any) => void;
  loading: boolean;
}

function DesignForm({ onSubmit, loading }: DesignFormProps) {
  const [prompt, setPrompt] = useState('');
  const [designType, setDesignType] = useState('analogue');
  const [referenceName, setReferenceName] = useState('');
  const [referenceSmiles, setReferenceSmiles] = useState('');
  const [targetName, setTargetName] = useState('');
  const [constraints, setConstraints] = useState({
    maxMw: 600,
    bbbPenetration: false,
    oralBioavailability: true,
  });

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    onSubmit({
      prompt,
      design_type: designType,
      reference_name: referenceName || undefined,
      reference_smiles: referenceSmiles || undefined,
      target_name: targetName || undefined,
      constraints: {
        max_mw: constraints.maxMw,
        bbb_penetration: constraints.bbbPenetration,
        oral_bioavailability: constraints.oralBioavailability,
      },
    });
  };

  return (
    <form onSubmit={handleSubmit} className="space-y-4">
      {/* Natural Language Prompt */}
      <div>
        <label className="block text-sm font-medium text-gray-700 mb-1">
          Design Request
        </label>
        <textarea
          value={prompt}
          onChange={(e) => setPrompt(e.target.value)}
          rows={3}
          className="w-full border rounded-md px-3 py-2 focus:ring-indigo-500 focus:border-indigo-500"
          placeholder="Describe the molecule you want to design..."
          required
        />
      </div>

      {/* Design Type */}
      <div>
        <label className="block text-sm font-medium text-gray-700 mb-1">
          Design Type
        </label>
        <select
          value={designType}
          onChange={(e) => setDesignType(e.target.value)}
          className="w-full border rounded-md px-3 py-2 focus:ring-indigo-500 focus:border-indigo-500"
        >
          <option value="analogue">Novel Analogue</option>
          <option value="bivalent">Bivalent Ligand</option>
          <option value="scaffold_hop">Scaffold Hop</option>
          <option value="de_novo">De Novo Design</option>
          <option value="optimization">Lead Optimization</option>
        </select>
      </div>

      {/* Reference Compound */}
      <div>
        <label className="block text-sm font-medium text-gray-700 mb-1">
          Reference Compound (optional)
        </label>
        <input
          type="text"
          value={referenceName}
          onChange={(e) => setReferenceName(e.target.value)}
          className="w-full border rounded-md px-3 py-2 focus:ring-indigo-500 focus:border-indigo-500"
          placeholder="e.g., KK103, Aspirin, etc."
        />
      </div>

      {/* Reference SMILES */}
      <div>
        <label className="block text-sm font-medium text-gray-700 mb-1">
          Reference SMILES (optional)
        </label>
        <input
          type="text"
          value={referenceSmiles}
          onChange={(e) => setReferenceSmiles(e.target.value)}
          className="w-full border rounded-md px-3 py-2 font-mono text-sm focus:ring-indigo-500 focus:border-indigo-500"
          placeholder="CC(=O)Oc1ccccc1C(=O)O"
        />
      </div>

      {/* Target Name */}
      <div>
        <label className="block text-sm font-medium text-gray-700 mb-1">
          Target (optional)
        </label>
        <input
          type="text"
          value={targetName}
          onChange={(e) => setTargetName(e.target.value)}
          className="w-full border rounded-md px-3 py-2 focus:ring-indigo-500 focus:border-indigo-500"
          placeholder="e.g., 5HT1A, mu-opioid receptor"
        />
      </div>

      {/* Constraints */}
      <div className="border-t pt-4">
        <h4 className="text-sm font-medium text-gray-700 mb-2">Constraints</h4>

        <div className="space-y-2">
          <div className="flex items-center justify-between">
            <label className="text-sm text-gray-600">Max MW</label>
            <input
              type="number"
              value={constraints.maxMw}
              onChange={(e) => setConstraints({
                ...constraints,
                maxMw: parseInt(e.target.value)
              })}
              className="w-24 border rounded px-2 py-1 text-sm"
            />
          </div>

          <div className="flex items-center justify-between">
            <label className="text-sm text-gray-600">BBB Penetration</label>
            <input
              type="checkbox"
              checked={constraints.bbbPenetration}
              onChange={(e) => setConstraints({
                ...constraints,
                bbbPenetration: e.target.checked
              })}
              className="rounded text-indigo-600"
            />
          </div>

          <div className="flex items-center justify-between">
            <label className="text-sm text-gray-600">Oral Bioavailability</label>
            <input
              type="checkbox"
              checked={constraints.oralBioavailability}
              onChange={(e) => setConstraints({
                ...constraints,
                oralBioavailability: e.target.checked
              })}
              className="rounded text-indigo-600"
            />
          </div>
        </div>
      </div>

      {/* Submit Button */}
      <button
        type="submit"
        disabled={loading || !prompt}
        className="w-full bg-indigo-600 text-white py-2 px-4 rounded-md hover:bg-indigo-700 disabled:opacity-50 disabled:cursor-not-allowed font-medium"
      >
        {loading ? 'Designing...' : 'Design Molecules'}
      </button>
    </form>
  );
}

export default DesignForm;
