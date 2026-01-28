import React, { useState } from 'react';
import DesignForm from './components/DesignForm';
import ResultsPanel from './components/ResultsPanel';
import MoleculeViewer from './components/MoleculeViewer';

interface DesignResult {
  job_id: string;
  status: string;
  molecules?: any[];
}

function App() {
  const [result, setResult] = useState<DesignResult | null>(null);
  const [selectedMolecule, setSelectedMolecule] = useState<any | null>(null);
  const [loading, setLoading] = useState(false);

  const handleDesignSubmit = async (formData: any) => {
    setLoading(true);
    try {
      const response = await fetch('/api/v1/design', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(formData),
      });
      const data = await response.json();

      // Poll for results
      const jobId = data.job_id;
      let jobResult = data;

      while (jobResult.status !== 'completed' && jobResult.status !== 'failed') {
        await new Promise(resolve => setTimeout(resolve, 2000));
        const statusResponse = await fetch(`/api/v1/design/${jobId}`);
        jobResult = await statusResponse.json();
      }

      setResult(jobResult);
      if (jobResult.molecules?.length > 0) {
        setSelectedMolecule(jobResult.molecules[0]);
      }
    } catch (error) {
      console.error('Design submission failed:', error);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-gray-100">
      {/* Header */}
      <header className="bg-indigo-600 text-white py-4 px-6 shadow-lg">
        <div className="max-w-7xl mx-auto flex justify-between items-center">
          <h1 className="text-2xl font-bold">Drug Design Platform</h1>
          <nav className="space-x-4">
            <a href="/docs" className="hover:text-indigo-200">API Docs</a>
          </nav>
        </div>
      </header>

      {/* Main Content */}
      <main className="max-w-7xl mx-auto py-6 px-4">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Design Form */}
          <div className="lg:col-span-1">
            <div className="bg-white rounded-lg shadow p-6">
              <h2 className="text-xl font-semibold mb-4">Design Request</h2>
              <DesignForm onSubmit={handleDesignSubmit} loading={loading} />
            </div>

            {/* Quick Actions */}
            <div className="bg-white rounded-lg shadow p-6 mt-6">
              <h3 className="text-lg font-semibold mb-3">Example Requests</h3>
              <div className="space-y-2">
                <button
                  onClick={() => handleDesignSubmit({
                    prompt: "Design a novel, patentable analogue of KK103 with improved metabolic stability",
                    reference_name: "KK103",
                    design_type: "analogue"
                  })}
                  className="w-full text-left px-3 py-2 bg-gray-50 rounded hover:bg-gray-100 text-sm"
                >
                  Novel KK103 analogue
                </button>
                <button
                  onClick={() => handleDesignSubmit({
                    prompt: "Design a novel dual 5HT1A-5HT2A agonist bivalent ligand",
                    design_type: "bivalent",
                    target_name: "5HT1A-5HT2A"
                  })}
                  className="w-full text-left px-3 py-2 bg-gray-50 rounded hover:bg-gray-100 text-sm"
                >
                  5HT1A-5HT2A bivalent ligand
                </button>
              </div>
            </div>
          </div>

          {/* Results Panel */}
          <div className="lg:col-span-2">
            {loading && (
              <div className="bg-white rounded-lg shadow p-6 text-center">
                <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-indigo-600 mx-auto"></div>
                <p className="mt-4 text-gray-600">Designing molecules...</p>
              </div>
            )}

            {result && !loading && (
              <div className="space-y-6">
                {/* Results Summary */}
                <div className="bg-white rounded-lg shadow p-6">
                  <h2 className="text-xl font-semibold mb-4">
                    Design Results
                    <span className="ml-2 text-sm font-normal text-gray-500">
                      {result.molecules?.length || 0} molecules designed
                    </span>
                  </h2>
                  <ResultsPanel
                    molecules={result.molecules || []}
                    onSelectMolecule={setSelectedMolecule}
                    selectedId={selectedMolecule?.smiles}
                  />
                </div>

                {/* Molecule Details */}
                {selectedMolecule && (
                  <div className="bg-white rounded-lg shadow p-6">
                    <h3 className="text-lg font-semibold mb-4">
                      Molecule Details - Rank #{selectedMolecule.rank}
                    </h3>
                    <MoleculeViewer molecule={selectedMolecule} />
                  </div>
                )}
              </div>
            )}

            {!result && !loading && (
              <div className="bg-white rounded-lg shadow p-6 text-center text-gray-500">
                <svg className="mx-auto h-12 w-12 text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
                </svg>
                <p className="mt-2">Submit a design request to get started</p>
              </div>
            )}
          </div>
        </div>
      </main>

      {/* Footer */}
      <footer className="bg-gray-800 text-gray-300 py-4 mt-8">
        <div className="max-w-7xl mx-auto px-4 text-center text-sm">
          Drug Design Platform - Automated Novel Molecule Design
        </div>
      </footer>
    </div>
  );
}

export default App;
