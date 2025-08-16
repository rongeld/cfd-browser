"use client";
import React, { useState } from "react";
import {
  ChevronDown,
  ChevronRight,
  BookOpen,
  Calculator,
  Eye,
  Lightbulb,
  Cpu,
} from "lucide-react";
import PDEStepVisualizer from "./PDEStepVisualizer";
import "katex/dist/katex.min.css";

// Simple inline math component to avoid external dependencies
const MathEquation = ({ children }: { children: string }) => {
  return (
    <span
      className="inline-block bg-gray-100 dark:bg-gray-700 px-2 py-1 rounded text-sm font-mono"
      style={{ fontFamily: "KaTeX_Math, Times New Roman, serif" }}
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
        style={{ fontFamily: "KaTeX_Math, Times New Roman, serif" }}
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

const CollapsibleSection = ({
  title,
  icon,
  children,
  defaultOpen = false,
}: CollapsibleSectionProps) => {
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
        {isOpen ? (
          <ChevronDown className="w-5 h-5" />
        ) : (
          <ChevronRight className="w-5 h-5" />
        )}
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
  simulationDebugData?: any;
  gridWidth?: number;
  gridHeight?: number;
}

const InteractiveVisualizer = ({
  onParameterChange,
}: InteractiveVisualizerProps) => {
  const [reynolds, setReynolds] = useState(33333); // Realistic starting value
  const [viscosity, setViscosity] = useState(0.000015); // Air viscosity
  const [velocity, setVelocity] = useState(10); // 10 m/s wind

  const handleReynoldsChange = (value: number) => {
    setReynolds(value);
    // Calculate corresponding viscosity: Re = vL/ν, so ν = vL/Re
    const characteristicLength = 0.5; // 50cm obstacle in 2m domain
    const newViscosity = (velocity * characteristicLength) / value;
    setViscosity(Math.max(0.000001, Math.min(0.001, newViscosity)));
    onParameterChange?.("viscosity", newViscosity);
  };

  const handleVelocityChange = (value: number) => {
    setVelocity(value);
    onParameterChange?.("windSpeed", value);
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
                  onParameterChange?.("viscosity", val);
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

const EducationalContent = ({
  onParameterChange,
  simulationDebugData,
  gridWidth = 20,
  gridHeight = 20,
}: InteractiveVisualizerProps) => {
  return (
    <div className="space-y-6">
      <div className="text-center mb-8">
        <h2 className="text-2xl font-bold mb-2">
          Learn CFD: Interactive Fluid Dynamics
        </h2>
        <p className="text-gray-600 dark:text-gray-400">
          Understand the physics and mathematics behind computational fluid
          dynamics
        </p>
      </div>

      <CollapsibleSection
        title="What is Computational Fluid Dynamics?"
        icon={<BookOpen className="w-5 h-5 text-blue-500" />}
        defaultOpen={true}
      >
        <div className="space-y-4">
          <p>
            Computational Fluid Dynamics (CFD) is a branch of fluid mechanics
            that uses numerical analysis to solve and analyze problems involving
            fluid flows. This simulator demonstrates how fluids behave when
            encountering obstacles.
          </p>

          <div className="bg-green-50 dark:bg-green-900/20 p-4 rounded-lg">
            <h4 className="font-semibold mb-2">Key Concepts:</h4>
            <ul className="space-y-2 text-sm">
              <li>
                <strong>Viscosity:</strong> A fluid's resistance to flow (like
                honey vs water)
              </li>
              <li>
                <strong>Pressure:</strong> Force per unit area exerted by the
                fluid
              </li>
              <li>
                <strong>Velocity Field:</strong> Speed and direction of fluid at
                every point
              </li>
              <li>
                <strong>Vorticity:</strong> Measure of fluid rotation and
                turbulence
              </li>
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
            The Navier-Stokes equations describe the motion of viscous fluids.
            They're based on Newton's second law applied to fluid motion.{" "}
            <strong>
              However, this simulator uses a simplified version for educational
              purposes!
            </strong>
          </p>

          <div className="space-y-4">
            <div>
              <h4 className="font-semibold mb-2">
                Full Navier-Stokes Momentum Equation:
              </h4>
              <BlockMath>∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f</BlockMath>
              <div className="text-sm text-gray-600 dark:text-gray-400">
                <p>
                  <strong>Where:</strong>
                </p>
                <ul className="mt-2 space-y-1 ml-4">
                  <li>
                    <MathEquation>∂u/∂t</MathEquation> = velocity change over
                    time
                  </li>
                  <li>
                    <MathEquation>(u·∇)u</MathEquation> = convection (nonlinear
                    term)
                  </li>
                  <li>
                    <MathEquation>∇p/ρ</MathEquation> = pressure gradient force
                  </li>
                  <li>
                    <MathEquation>ν∇²u</MathEquation> = viscous diffusion
                  </li>
                  <li>
                    <MathEquation>f</MathEquation> = external forces (gravity,
                    etc.)
                  </li>
                </ul>
              </div>
            </div>

            <div>
              <h4 className="font-semibold mb-2">
                Continuity Equation (Mass Conservation):
              </h4>
              <BlockMath>∇·u = 0</BlockMath>
              <p className="text-sm text-gray-600 dark:text-gray-400">
                This ensures that mass is conserved - what flows in must flow
                out (incompressible flow).
              </p>
            </div>

            <div className="bg-orange-50 dark:bg-orange-900/20 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">
                What This Simulator Actually Solves:
              </h4>
              <div className="text-sm space-y-2">
                <p>
                  <strong>1. Operator Splitting:</strong> Instead of solving
                  everything at once, we split into steps:
                </p>
                <ul className="ml-4 space-y-1">
                  <li>
                    • <strong>Advection step:</strong> Solve ∂u/∂t + (u·∇)u = 0
                  </li>
                  <li>
                    • <strong>Viscosity step:</strong> Solve ∂u/∂t = ν∇²u
                  </li>
                  <li>
                    • <strong>Projection step:</strong> Solve ∇²p = -∇·u, then u
                    = u - ∇p
                  </li>
                </ul>
                <p>
                  <strong>2. Simplified assumptions:</strong> No body forces,
                  constant properties, 2D only
                </p>
              </div>
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
            The Reynolds number determines the flow behavior - whether it's
            smooth (laminar) or chaotic (turbulent). You're absolutely right
            that Reynolds numbers are often in the thousands!
          </p>

          <BlockMath>Re = ρvL/μ = vL/ν</BlockMath>

          <div className="text-sm text-gray-600 dark:text-gray-400 mb-4">
            <ul className="space-y-1 ml-4">
              <li>
                <MathEquation>ρ</MathEquation> = fluid density
              </li>
              <li>
                <MathEquation>v</MathEquation> = characteristic velocity
              </li>
              <li>
                <MathEquation>L</MathEquation> = characteristic length
              </li>
              <li>
                <MathEquation>μ</MathEquation> = dynamic viscosity
              </li>
            </ul>
          </div>

          <div className="bg-yellow-50 dark:bg-yellow-900/20 p-4 rounded-lg mb-4">
            <h4 className="font-semibold mb-2">Real-World Examples:</h4>
            <ul className="text-sm space-y-1">
              <li>
                <strong>Re ≈ 0.1:</strong> Bacteria swimming in water
              </li>
              <li>
                <strong>Re ≈ 100:</strong> Small insects flying
              </li>
              <li>
                <strong>Re ≈ 1,000:</strong> Fish swimming
              </li>
              <li>
                <strong>Re ≈ 10,000:</strong> Cyclist in air
              </li>
              <li>
                <strong>Re ≈ 100,000:</strong> Car moving at highway speed
              </li>
              <li>
                <strong>Re ≈ 1,000,000:</strong> Commercial aircraft
              </li>
              <li>
                <strong>Re ≈ 10,000,000+:</strong> Large ships
              </li>
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
          <div className="bg-red-50 dark:bg-red-900/20 p-4 rounded-lg">
            <h4 className="font-semibold mb-2">
              ⚠️ Important: This is NOT Full Navier-Stokes
            </h4>
            <p className="text-sm">
              This simulator uses a{" "}
              <strong>simplified incompressible flow solver</strong> - not the
              complete Navier-Stokes equations. It's perfect for learning basic
              concepts but missing many real-world complexities.
            </p>
          </div>

          <div className="space-y-4">
            <h4 className="font-semibold">What's Actually Implemented:</h4>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="bg-blue-50 dark:bg-blue-800 p-4 rounded-lg">
                <h5 className="font-semibold mb-2">
                  ✅ Pressure Projection Method
                </h5>
                <p className="text-sm">
                  Splits each timestep into: advection → viscosity → pressure
                  projection. Uses Jacobi iterations to solve the pressure
                  Poisson equation.
                </p>
              </div>

              <div className="bg-green-50 dark:bg-green-800 p-4 rounded-lg">
                <h5 className="font-semibold mb-2">
                  ✅ Semi-Lagrangian Advection
                </h5>
                <p className="text-sm">
                  Handles the nonlinear term (u·∇)u by tracing particles
                  backward in time. Stable but diffusive.
                </p>
              </div>

              <div className="bg-purple-50 dark:bg-purple-800 p-4 rounded-lg">
                <h5 className="font-semibold mb-2">✅ Explicit Viscosity</h5>
                <p className="text-sm">
                  Simple diffusion using finite differences. Fast but can be
                  unstable with high viscosity.
                </p>
              </div>

              <div className="bg-yellow-50 dark:bg-yellow-800 p-4 rounded-lg">
                <h5 className="font-semibold mb-2">⚠️ Simplified Boundaries</h5>
                <p className="text-sm">
                  Basic no-slip conditions on obstacles. Real CFD uses more
                  sophisticated boundary treatments.
                </p>
              </div>
            </div>

            <div className="bg-gray-50 dark:bg-gray-800 p-4 rounded-lg">
              <h5 className="font-semibold mb-2">
                Missing from Full Navier-Stokes:
              </h5>
              <ul className="text-sm space-y-1 list-disc ml-4">
                <li>Body forces (gravity, magnetic fields)</li>
                <li>Variable density and viscosity</li>
                <li>Turbulence modeling</li>
                <li>Heat transfer coupling</li>
                <li>Compressibility effects</li>
                <li>Advanced boundary conditions</li>
                <li>High-order spatial discretizations</li>
              </ul>
            </div>
          </div>
        </div>
      </CollapsibleSection>

      <CollapsibleSection
        title="Live PDE Solution Process"
        icon={<Cpu className="w-5 h-5 text-cyan-500" />}
      >
        <div className="space-y-4">
          <p>
            Watch how the computer solves the fluid equations step by step! This
            shows the actual numerical process happening behind the scenes.
          </p>

          <PDEStepVisualizer
            gridWidth={gridWidth}
            gridHeight={gridHeight}
            simulationData={simulationDebugData}
            onStepUpdate={(step, data) => {
              console.log(`Step ${step} updated:`, data);
            }}
          />
        </div>
      </CollapsibleSection>

      <CollapsibleSection
        title="Real CFD vs This Simulator"
        icon={<Calculator className="w-5 h-5 text-red-500" />}
      >
        <div className="space-y-4">
          <p>
            This simulator is great for learning, but real CFD software is much
            more sophisticated! Here's what professional tools like ANSYS
            Fluent, OpenFOAM, or SU2 actually do:
          </p>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="bg-green-50 dark:bg-green-900/20 p-4 rounded-lg">
              <h4 className="font-semibold mb-2 text-green-700 dark:text-green-300">
                🎓 This Educational Simulator
              </h4>
              <ul className="text-sm space-y-1">
                <li>• 2D rectangular grid (80-300 cells)</li>
                <li>• Operator splitting method</li>
                <li>• Explicit time stepping</li>
                <li>• Simple boundary conditions</li>
                <li>• No turbulence modeling</li>
                <li>• Fixed fluid properties</li>
                <li>• ~60 FPS real-time</li>
              </ul>
            </div>

            <div className="bg-blue-50 dark:bg-blue-900/20 p-4 rounded-lg">
              <h4 className="font-semibold mb-2 text-blue-700 dark:text-blue-300">
                🏭 Professional CFD Software
              </h4>
              <ul className="text-sm space-y-1">
                <li>• 3D unstructured meshes (millions+ cells)</li>
                <li>• Coupled/segregated solvers</li>
                <li>• Implicit time stepping</li>
                <li>• Complex boundary conditions</li>
                <li>• Advanced turbulence models (RANS, LES, DNS)</li>
                <li>• Variable properties, heat transfer</li>
                <li>• Hours to days for convergence</li>
              </ul>
            </div>
          </div>

          <div className="bg-yellow-50 dark:bg-yellow-900/20 p-4 rounded-lg">
            <h4 className="font-semibold mb-2">Real CFD Equation Systems:</h4>
            <div className="text-sm space-y-2">
              <p>
                <strong>RANS (Reynolds-Averaged Navier-Stokes):</strong>
              </p>
              <div className="bg-white dark:bg-gray-800 p-2 rounded font-mono text-xs">
                ∂(ρū)/∂t + ∇·(ρū⊗ū) = -∇p̄ + ∇·(μ∇ū) + ∇·τ_turb + f
              </div>
              <p className="text-xs">
                Where τ_turb requires additional turbulence equations (k-ε, k-ω,
                SST, etc.)
              </p>

              <p className="mt-2">
                <strong>Plus energy equation for heat transfer:</strong>
              </p>
              <div className="bg-white dark:bg-gray-800 p-2 rounded font-mono text-xs">
                ∂(ρcₚT)/∂t + ∇·(ρcₚūT) = ∇·(λ∇T) + Φ
              </div>
            </div>
          </div>

          <div className="bg-gray-50 dark:bg-gray-800 p-4 rounded-lg">
            <h4 className="font-semibold mb-2">
              Why This Matters for Learning:
            </h4>
            <p className="text-sm">
              Understanding this simplified version gives you the foundation to
              understand real CFD! The basic concepts (discretization,
              pressure-velocity coupling, boundary conditions) are the same -
              just much more sophisticated in professional tools.
            </p>
          </div>
        </div>
      </CollapsibleSection>

      <CollapsibleSection
        title="Experiment Ideas"
        icon={<Lightbulb className="w-5 h-5 text-yellow-500" />}
      >
        <div className="space-y-4">
          <p>Try these experiments to understand different fluid phenomena:</p>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="border border-gray-200 dark:border-gray-700 p-4 rounded-lg">
              <h4 className="font-semibold mb-2">🔵 Cylinder Flow</h4>
              <p className="text-sm mb-2">
                Draw a circular obstacle and observe:
              </p>
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
