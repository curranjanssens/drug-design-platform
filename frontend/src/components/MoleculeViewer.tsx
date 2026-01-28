import React, { useEffect, useRef } from 'react';

interface MoleculeViewerProps {
  molecule: any;
}

function MoleculeViewer({ molecule }: MoleculeViewerProps) {
  const viewerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    // 3Dmol.js integration would go here
    // For now, we'll just display the structure data
  }, [molecule]);

  return (
    <div className="space-y-6">
      {/* Structure Display */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        {/* 2D Structure */}
        <div>
          <h4 className="text-sm font-medium text-gray-700 mb-2">2D Structure</h4>
          <div className="bg-gray-50 rounded-lg p-4 h-48 flex items-center justify-center">
            {molecule.image_svg ? (
              <div dangerouslySetInnerHTML={{ __html: molecule.image_svg }} />
            ) : (
              <div className="text-gray-400 text-sm">
                Structure image not available
              </div>
            )}
          </div>
        </div>

        {/* 3D Viewer Placeholder */}
        <div>
          <h4 className="text-sm font-medium text-gray-700 mb-2">3D Structure</h4>
          <div
            ref={viewerRef}
            className="bg-gray-900 rounded-lg h-48 flex items-center justify-center"
          >
            <div className="text-gray-400 text-sm">
              3D viewer (requires 3Dmol.js)
            </div>
          </div>
        </div>
      </div>

      {/* SMILES */}
      <div>
        <h4 className="text-sm font-medium text-gray-700 mb-1">SMILES</h4>
        <div className="bg-gray-50 rounded px-3 py-2 font-mono text-sm break-all">
          {molecule.smiles}
        </div>
      </div>

      {/* Key Properties */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        <PropertyCard
          label="Molecular Weight"
          value={molecule.properties?.molecular_weight?.toFixed(1)}
          unit="Da"
        />
        <PropertyCard
          label="LogP"
          value={molecule.properties?.logp?.toFixed(2)}
        />
        <PropertyCard
          label="QED"
          value={molecule.properties?.qed?.toFixed(3)}
        />
        <PropertyCard
          label="TPSA"
          value={molecule.properties?.tpsa?.toFixed(1)}
          unit="A2"
        />
      </div>

      {/* Novelty Assessment */}
      <div>
        <h4 className="text-sm font-medium text-gray-700 mb-2">Novelty Assessment</h4>
        <div className="bg-gray-50 rounded-lg p-4 grid grid-cols-2 gap-4">
          <div>
            <span className="text-xs text-gray-500">Similarity to Nearest</span>
            <p className="font-medium">
              {((molecule.novelty?.tanimoto_similarity || 0) * 100).toFixed(1)}%
            </p>
          </div>
          <div>
            <span className="text-xs text-gray-500">Patentable</span>
            <p className={`font-medium ${
              molecule.novelty?.patentable ? 'text-green-600' : 'text-red-600'
            }`}>
              {molecule.novelty?.patentable ? 'Yes' : 'No'}
            </p>
          </div>
          <div>
            <span className="text-xs text-gray-500">FTO Risk</span>
            <p className={`font-medium ${
              molecule.novelty?.fto_risk === 'low' ? 'text-green-600' :
              molecule.novelty?.fto_risk === 'medium' ? 'text-yellow-600' : 'text-red-600'
            }`}>
              {molecule.novelty?.fto_risk?.toUpperCase() || 'Unknown'}
            </p>
          </div>
          <div>
            <span className="text-xs text-gray-500">Scaffold Novel</span>
            <p className={`font-medium ${
              molecule.novelty?.scaffold_novel ? 'text-green-600' : 'text-yellow-600'
            }`}>
              {molecule.novelty?.scaffold_novel ? 'Yes' : 'No'}
            </p>
          </div>
        </div>
      </div>

      {/* ADMET Summary */}
      <div>
        <h4 className="text-sm font-medium text-gray-700 mb-2">ADMET Profile</h4>
        <div className="bg-gray-50 rounded-lg p-4">
          <div className="flex items-center justify-between mb-2">
            <span className="text-sm">Safety Score</span>
            <span className={`font-medium ${
              (molecule.admet_safety_score || 0) >= 0.7 ? 'text-green-600' :
              (molecule.admet_safety_score || 0) >= 0.5 ? 'text-yellow-600' : 'text-red-600'
            }`}>
              {((molecule.admet_safety_score || 0) * 100).toFixed(0)}%
            </span>
          </div>
          <div className="w-full bg-gray-200 rounded-full h-2">
            <div
              className={`h-2 rounded-full ${
                (molecule.admet_safety_score || 0) >= 0.7 ? 'bg-green-500' :
                (molecule.admet_safety_score || 0) >= 0.5 ? 'bg-yellow-500' : 'bg-red-500'
              }`}
              style={{ width: `${(molecule.admet_safety_score || 0) * 100}%` }}
            />
          </div>
          {molecule.admet?.warnings?.length > 0 && (
            <div className="mt-3">
              <span className="text-xs text-gray-500">Warnings</span>
              <ul className="mt-1 text-sm text-red-600">
                {molecule.admet.warnings.map((w: string, i: number) => (
                  <li key={i}>• {w}</li>
                ))}
              </ul>
            </div>
          )}
        </div>
      </div>

      {/* Synthesis Route */}
      {molecule.synthesis_steps && (
        <div>
          <h4 className="text-sm font-medium text-gray-700 mb-2">Synthesis Overview</h4>
          <div className="bg-gray-50 rounded-lg p-4">
            <div className="grid grid-cols-3 gap-4 text-center">
              <div>
                <span className="text-2xl font-bold text-indigo-600">
                  {molecule.synthesis_steps}
                </span>
                <p className="text-xs text-gray-500">Steps</p>
              </div>
              <div>
                <span className="text-2xl font-bold text-indigo-600">
                  {molecule.best_route?.overall_yield
                    ? `${(molecule.best_route.overall_yield * 100).toFixed(0)}%`
                    : 'N/A'}
                </span>
                <p className="text-xs text-gray-500">Est. Yield</p>
              </div>
              <div>
                <span className="text-2xl font-bold text-indigo-600">
                  {molecule.best_route?.confidence
                    ? `${(molecule.best_route.confidence * 100).toFixed(0)}%`
                    : 'N/A'}
                </span>
                <p className="text-xs text-gray-500">Confidence</p>
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

function PropertyCard({
  label,
  value,
  unit,
}: {
  label: string;
  value: string | number | undefined;
  unit?: string;
}) {
  return (
    <div className="bg-gray-50 rounded-lg p-3 text-center">
      <p className="text-xs text-gray-500">{label}</p>
      <p className="text-lg font-semibold text-gray-900">
        {value || 'N/A'}
        {unit && <span className="text-xs text-gray-500 ml-1">{unit}</span>}
      </p>
    </div>
  );
}

export default MoleculeViewer;
