import React from 'react';

interface Molecule {
  rank: number;
  smiles: string;
  inchi_key: string;
  score: number;
  properties: {
    molecular_weight: number;
    logp: number;
    qed: number;
  };
  novelty: {
    tanimoto_similarity: number;
    patentable: boolean;
    fto_risk: string;
  };
  admet_safety_score: number;
  synthesis_steps: number;
}

interface ResultsPanelProps {
  molecules: Molecule[];
  onSelectMolecule: (molecule: Molecule) => void;
  selectedId: string | null;
}

function ResultsPanel({ molecules, onSelectMolecule, selectedId }: ResultsPanelProps) {
  if (molecules.length === 0) {
    return (
      <div className="text-center text-gray-500 py-8">
        No molecules generated yet
      </div>
    );
  }

  return (
    <div className="overflow-x-auto">
      <table className="min-w-full divide-y divide-gray-200">
        <thead className="bg-gray-50">
          <tr>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              Rank
            </th>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              Score
            </th>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              MW
            </th>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              QED
            </th>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              Novelty
            </th>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              Patentable
            </th>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              Safety
            </th>
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase">
              Steps
            </th>
          </tr>
        </thead>
        <tbody className="bg-white divide-y divide-gray-200">
          {molecules.map((mol) => (
            <tr
              key={mol.smiles}
              onClick={() => onSelectMolecule(mol)}
              className={`cursor-pointer hover:bg-gray-50 ${
                selectedId === mol.smiles ? 'bg-indigo-50' : ''
              }`}
            >
              <td className="px-3 py-2 whitespace-nowrap text-sm font-medium text-gray-900">
                #{mol.rank}
              </td>
              <td className="px-3 py-2 whitespace-nowrap text-sm">
                <span className={`font-medium ${
                  mol.score >= 0.7 ? 'text-green-600' :
                  mol.score >= 0.5 ? 'text-yellow-600' : 'text-red-600'
                }`}>
                  {(mol.score * 100).toFixed(0)}%
                </span>
              </td>
              <td className="px-3 py-2 whitespace-nowrap text-sm text-gray-500">
                {mol.properties?.molecular_weight?.toFixed(0) || 'N/A'}
              </td>
              <td className="px-3 py-2 whitespace-nowrap text-sm text-gray-500">
                {mol.properties?.qed?.toFixed(2) || 'N/A'}
              </td>
              <td className="px-3 py-2 whitespace-nowrap text-sm text-gray-500">
                {((1 - (mol.novelty?.tanimoto_similarity || 0)) * 100).toFixed(0)}%
              </td>
              <td className="px-3 py-2 whitespace-nowrap text-sm">
                {mol.novelty?.patentable ? (
                  <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium bg-green-100 text-green-800">
                    Yes
                  </span>
                ) : (
                  <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium bg-red-100 text-red-800">
                    No
                  </span>
                )}
              </td>
              <td className="px-3 py-2 whitespace-nowrap text-sm">
                <span className={`font-medium ${
                  (mol.admet_safety_score || 0) >= 0.7 ? 'text-green-600' :
                  (mol.admet_safety_score || 0) >= 0.5 ? 'text-yellow-600' : 'text-red-600'
                }`}>
                  {((mol.admet_safety_score || 0) * 100).toFixed(0)}%
                </span>
              </td>
              <td className="px-3 py-2 whitespace-nowrap text-sm text-gray-500">
                {mol.synthesis_steps || 'N/A'}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

export default ResultsPanel;
