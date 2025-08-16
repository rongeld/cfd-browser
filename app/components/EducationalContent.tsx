"use client";
import React, { useState } from 'react';
import { ChevronDown, ChevronRight, BookOpen, Calculator, Eye, Lightbulb } from 'lucide-react';
import 'katex/dist/katex.min.css';

// Simple inline math component to avoid external dependencies
const MathEquation = ({ children }: { children: string }) => {
  return (
    <span 
      className="inline-block bg-gray-100 dark:bg-gray-700 px-2 py-1 rounded text-sm font-mono"
      style={{ fontFamily: 'KaTeX_Math, Times New Roman, serif' }}
    >
      {children}
    </span>
  );
};

const BlockMath = ({ children }: { children: string }) => {
  return (
    <div className="bg-gray-100 dark:bg-gray-700 p-4 rounded-lg my-4 text-center">
      <span 
        className="text-lg font-mono"
        style={{ fontFamily: 'KaTeX_Math, Times New Roman, serif' }}
      >
        {children}
      </span>
    </div>
  );
};

interface CollapsibleSectionProps {
  title: string;
  icon: React.ReactNode;
  children: React.ReactNode;
  defaultOpen?: boolean;
}

const CollapsibleSection = ({ title, icon, children, defaultOpen = false }: CollapsibleSectionProps) => {
  const [isOpen, setIsOpen] = useState(defaultOpen);

  return (
    <div className="border border-gray-200 dark:border-gray-700 rounded-lg mb-4">
      <button
        onClick={() => setIsOpen(!isOpen)}
        className="w-full flex items-center justify-between p-4 text-left hover:bg-gray-50 dark:hover:bg-gray-800 transition-colors"
      >
        <div className="flex items-center gap-3">
          {icon}
          <h3 className="text-lg font-semibold">{title}</h3>
        </div>
        {isOpen ? <ChevronDown className="w-5 h-5" /> : <ChevronRight className="w-5 h-5" />}
      </button>
      {isOpen && (
        <div className="px-4 pb-4 border-t border-gray-200 dark:border-gray-700">
          {children}
        </div>
      )}
    </div>
  );
};

interface InteractiveVisualizerProps {
  onParameterChange?: (param: string, value: number) => void;
}

const InteractiveVisualizer = ({ onParameterChange }: InteractiveVisualizerProps) => {
  const [reynolds, setReynolds] = useState(33333); // Realistic starting value
  const [viscosity, setViscosity] = useState(0.000015); // Air viscosity
  const [velocity, setVelocity] = useState(10); // 10 m/s wind

  const handleReynoldsChange = (value: number) => {
    setReynolds(value);
    // Calculate corresponding viscosity: Re = vL/ν, so ν = vL/Re
    const characteristicLength = 0.5; // 50cm obstacle in 2m domain
    const newViscosity = (velocity * characteristicLength) / value;
    setViscosity(Math.max(0.000001, Math.min(0.001, newViscosity)));
    onParameterChange?.('viscosity', newViscosity);
  };

  const handleVelocityChange = (value: number) => {
    setVelocity(value);
    onParameterChange?.('windSpeed', value);
    // Recalculate Reynolds with new velocity
    const characteristicLength = 0.5;
    const newReynolds = (value * characteristicLength) / viscosity;
    setReynolds(newReynolds);
  };

  const getFlowDescription = (re: number) => {
    if (re < 1) return "Creeping flow - viscous forces dominate";
    if (re < 40) return "Steady laminar flow - predictable streamlines";
    if (re < 200) return "Vortex shedding begins - oscillatory behavior";
    if (re < 1000) return "Turbulent wake - chaotic flow patterns";
    if (re < 10000) return "Fully turbulent - complex vortex structures";
    return "Highly turbulent - extreme flow complexity";
  };

  return (
    <div className="space-y-4">
      <div className="bg-blue-50 dark:bg-blue-900/20 p-4 rounded-lg">
        <h4 className="font-semibold mb-2">Reynolds Number Calculator</h4>
        <p className="text-sm text-gray-600 dark:text-gray-400 mb-4">
          Adjust parameters to see how Reynolds number affects flow behavior
        </p>
        
        <div className="space-y-4">
          <div>
            <label className="block text-sm font-medium mb-2">
              Reynolds Number: {reynolds.toFixed(0)}
            </label>
            <input
              type="range"
              min="1"
              max="50000"
              step="50"
              value={reynolds}
              onChange={(e) => handleReynoldsChange(parseInt(e.target.value))}
              className="w-full"
            />
            <p className="text-xs mt-1 text-blue-600 dark:text-blue-400">
              {getFlowDescription(reynolds)}
            </p>
          </div>
          
          <div className="grid grid-cols-2 gap-4">
            <div>
              <label className="block text-sm font-medium mb-1">
                Viscosity: {viscosity.toExponential(2)} m²/s
              </label>
              <input
                type="range"
                min="0.000001"
                max="0.0001"
                step="0.000001"
                value={viscosity}
                onChange={(e) => {
                  const val = parseFloat(e.target.value);
                  setViscosity(val);
                  onParameterChange?.('viscosity', val);
                  // Recalculate Reynolds
                  const characteristicLength = 0.5;
                  const newReynolds = (velocity * characteristicLength) / val;
                  setReynolds(newReynolds);
                }}
                className="w-full"
              />
              <p className="text-xs text-gray-500 mt-1">
                Air: 1.5×10⁻⁵ | Water: 1×10⁻⁶
              </p>
            </div>
            
            <div>
              <label className="block text-sm font-medium mb-1">
                Velocity: {velocity.toFixed(1)} m/s
              </label>
              <input
                type="range"
                min="1"
                max="50"
                value={velocity}
                onChange={(e) => handleVelocityChange(parseInt(e.target.value))}
                className="w-full"
              />
              <p className="text-xs text-gray-500 mt-1">
                Domain: 2m × 1m | Obstacle: ~0.5m
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

const EducationalContent = ({ onParameterChange }: InteractiveVisualizerProps) => {
  return (
    <div className="space-y-6">
      <div className="text-center mb-8">
        <h2 className="text-2xl font-bold mb-2">Learn CFD: Interactive Fluid Dynamics</h2>
        <p className="text-gray-600 dark:text-gray-400">
          Understand the physics and mathematics behind computational fluid dynamics
        </p>
      </div>

      <CollapsibleSection 
        title="What is Computational Fluid Dynamics?" 
        icon={<BookOpen className="w-5 h-5 text-blue-500" />}
        defaultOpen={true}
      >
        <div className="space-y-4">
          <p>
            Computational Fluid Dynamics (CFD) is a branch of fluid mechanics that uses numerical analysis 
            to solve and analyze problems involving fluid flows. This simulator demonstrates how fluids 
            behave when encountering obstacles.
          </p>
          
          <div className="bg-green-50 dark:bg-green-900/20 p-4 rounded-lg">
            <h4 className="font-semibold mb-2">Key Concepts:</h4>
            <ul className="space-y-2 text-sm">
              <li><strong>Viscosity:</strong> A fluid's resistance to flow (like honey vs water)</li>
              <li><strong>Pressure:</strong> Force per unit area exerted by the fluid</li>
              <li><strong>Velocity Field:</strong> Speed and direction of fluid at every point</li>
              <li><strong>Vorticity:</strong> Measure of fluid rotation and turbulence</li>
            </ul>
          </div>
        </div>
      </CollapsibleSection>

      <CollapsibleSection 
        title="The Navier-Stokes Equations" 
        icon={<Calculator className="w-5 h-5 text-green-500" />}
      >
        <div className="space-y-4">
          <p>
            The Navier-Stokes equations describe the motion of viscous fluids. They're based on 
            Newton's second law applied to fluid motion.
          </p>
          
          <div className="space-y-4">
            <div>
              <h4 className="font-semibold mb-2">Momentum Equation:</h4>
              <BlockMath>∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f</BlockMath>
              <div className="text-sm text-gray-600 dark:text-gray-400">
                <p><strong>Where:</strong></p>
                <ul className="mt-2 space-y-1 ml-4">
                  <li><MathEquation>u</MathEquation> = velocity vector</li>
                  <li><MathEquation>p</MathEquation> = pressure</li>
                  <li><MathEquation>ρ</MathEquation> = density</li>
                  <li><MathEquation>ν</MathEquation> = kinematic viscosity</li>
                  <li><MathEquation>f</MathEquation> = external forces</li>
                </ul>
              </div>
            </div>
            
            <div>
              <h4 className="font-semibold mb-2">Continuity Equation (Mass Conservation):</h4>
              <BlockMath>∇·u = 0</BlockMath>
              <p className="text-sm text-gray-600 dark:text-gray-400">
                This ensures that mass is conserved - what flows in must flow out.
              </p>
            </div>
          </div>
        </div>
      </CollapsibleSection>

      <CollapsibleSection 
        title="Reynolds Number & Flow Regimes" 
        icon={<Calculator className="w-5 h-5 text-purple-500" />}
      >
        <div className="space-y-4">
          <p>
            The Reynolds number determines the flow behavior - whether it's smooth (laminar) 
            or chaotic (turbulent). You're absolutely right that Reynolds numbers are often in the thousands!
          </p>
          
          <BlockMath>Re = ρvL/μ = vL/ν</BlockMath>
          
          <div className="text-sm text-gray-600 dark:text-gray-400 mb-4">
            <ul className="space-y-1 ml-4">
              <li><MathEquation>ρ</MathEquation> = fluid density</li>
              <li><MathEquation>v</MathEquation> = characteristic velocity</li>
              <li><MathEquation>L</MathEquation> = characteristic length</li>
              <li><MathEquation>μ</MathEquation> = dynamic viscosity</li>
            </ul>
          </div>

          <div className="bg-yellow-50 dark:bg-yellow-900/20 p-4 rounded-lg mb-4">
            <h4 className="font-semibold mb-2">Real-World Examples:</h4>
            <ul className="text-sm space-y-1">
              <li><strong>Re ≈ 0.1:</strong> Bacteria swimming in water</li>
              <li><strong>Re ≈ 100:</strong> Small insects flying</li>
              <li><strong>Re ≈ 1,000:</strong> Fish swimming</li>
              <li><strong>Re ≈ 10,000:</strong> Cyclist in air</li>
              <li><strong>Re ≈ 100,000:</strong> Car moving at highway speed</li>
              <li><strong>Re ≈ 1,000,000:</strong> Commercial aircraft</li>
              <li><strong>Re ≈ 10,000,000+:</strong> Large ships</li>
            </ul>
          </div>
          
          <InteractiveVisualizer onParameterChange={onParameterChange} />
        </div>
      </CollapsibleSection>

      <CollapsibleSection 
        title="Numerical Methods" 
        icon={<Eye className="w-5 h-5 text-orange-500" />}
      >
        <div className="space-y-4">
          <p>
            Since the Navier-Stokes equations are too complex to solve analytically, 
            we use numerical methods to approximate the solution.
          </p>
          
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="bg-gray-50 dark:bg-gray-800 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">Finite Difference Method</h4>
              <p className="text-sm">
                Approximates derivatives using discrete grid points. Simple but can be 
                less accurate for complex geometries.
              </p>
            </div>
            
            <div className="bg-gray-50 dark:bg-gray-800 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">Lattice Boltzmann Method</h4>
              <p className="text-sm">
                Models fluid as particles moving on a lattice. Great for complex boundaries 
                and parallel computation.
              </p>
            </div>
          </div>
          
          <div className="bg-blue-50 dark:bg-blue-900/20 p-4 rounded-lg">
            <h4 className="font-semibold mb-2">This Simulator Uses:</h4>
            <ul className="text-sm space-y-1">
              <li>• Finite difference discretization on a regular grid</li>
              <li>• Explicit time stepping for velocity updates</li>
              <li>• Pressure projection to enforce incompressibility</li>
              <li>• Particle tracing for flow visualization</li>
            </ul>
          </div>
        </div>
      </CollapsibleSection>

      <CollapsibleSection 
        title="Experiment Ideas" 
        icon={<Lightbulb className="w-5 h-5 text-yellow-500" />}
      >
        <div className="space-y-4">
          <p>
            Try these experiments to understand different fluid phenomena:
          </p>
          
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="border border-gray-200 dark:border-gray-700 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">🔵 Cylinder Flow</h4>
              <p className="text-sm mb-2">Draw a circular obstacle and observe:</p>
              <ul className="text-xs space-y-1">
                <li>• Low Re: Symmetric streamlines</li>
                <li>• Medium Re: Vortex shedding (Kármán street)</li>
                <li>• High Re: Turbulent wake</li>
              </ul>
            </div>
            
            <div className="border border-gray-200 dark:border-gray-700 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">🟪 Venturi Effect</h4>
              <p className="text-sm mb-2">Create a channel constriction:</p>
              <ul className="text-xs space-y-1">
                <li>• Velocity increases in narrow section</li>
                <li>• Pressure decreases (Bernoulli's principle)</li>
                <li>• Watch pressure field visualization</li>
              </ul>
            </div>
            
            <div className="border border-gray-200 dark:border-gray-700 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">🔺 Airfoil Shape</h4>
              <p className="text-sm mb-2">Draw wing-like shapes:</p>
              <ul className="text-xs space-y-1">
                <li>• Asymmetric flow creates lift</li>
                <li>• Higher velocity on top surface</li>
                <li>• Lower pressure above (Coanda effect)</li>
              </ul>
            </div>
            
            <div className="border border-gray-200 dark:border-gray-700 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">⚡ Multiple Obstacles</h4>
              <p className="text-sm mb-2">Place several objects:</p>
              <ul className="text-xs space-y-1">
                <li>• Interference patterns</li>
                <li>• Accelerated flow between objects</li>
                <li>• Complex vortex interactions</li>
              </ul>
            </div>
          </div>
        </div>
      </CollapsibleSection>
    </div>
  );
};

export default EducationalContent;
