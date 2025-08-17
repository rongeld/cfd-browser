"use client";
import React, { useState, useEffect, useRef } from "react";
import { ChevronDown, ChevronRight } from "lucide-react";
import katex from "katex";
import "katex/dist/katex.min.css";

// Custom KaTeX component
const Math: React.FC<{ children: string; block?: boolean }> = ({
  children,
  block = false,
}) => {
  const ref = useRef<HTMLSpanElement>(null);

  useEffect(() => {
    if (ref.current) {
      try {
        katex.render(children, ref.current, {
          displayMode: block,
          throwOnError: false,
        });
      } catch (error) {
        console.warn("KaTeX rendering error:", error);
        if (ref.current) {
          ref.current.textContent = children;
        }
      }
    }
  }, [children, block]);

  return <span ref={ref} className={block ? "block my-2" : "inline"} />;
};

interface FAQItem {
  question: string;
  answer: string | React.ReactNode;
  category: string;
}

interface FAQSectionProps {
  className?: string;
}

const FAQSection: React.FC<FAQSectionProps> = ({ className = "" }) => {
  const [openItems, setOpenItems] = useState<Set<number>>(new Set());
  const [activeCategory, setActiveCategory] = useState("all");

  const toggleItem = (index: number) => {
    const newOpenItems = new Set(openItems);
    if (newOpenItems.has(index)) {
      newOpenItems.delete(index);
    } else {
      newOpenItems.add(index);
    }
    setOpenItems(newOpenItems);
  };

  const faqData: FAQItem[] = [
    // Overview & Basics
    {
      category: "overview",
      question: "What is this CFD Simulator and what does it do?",
      answer: (
        <div>
          <p className="mb-2">
            This is a real-time, interactive CFD (Computational Fluid Dynamics)
            simulator built with React/Next.js that solves the incompressible
            Navier-Stokes equations using numerical methods. It features
            boundary layer visualization, particle tracing, pressure field
            analysis, and interactive obstacle drawing.
          </p>
          <p className="mb-2">
            The simulator implements a 2D incompressible Navier-Stokes solver
            using:
          </p>
          <ul className="list-disc pl-5 space-y-1">
            <li>
              <strong>Projection Method</strong> for pressure-velocity coupling
            </li>
            <li>
              <strong>Semi-Lagrangian Advection</strong> for stability at high
              Reynolds numbers
            </li>
            <li>
              <strong>Explicit Viscous Diffusion</strong> with boundary layer
              enhancement
            </li>
            <li>
              <strong>Real-time Visualization</strong> with multiple rendering
              modes
            </li>
          </ul>
        </div>
      ),
    },
    {
      category: "overview",
      question: "What are the key equations being solved?",
      answer: (
        <div>
          <p className="mb-2">
            The simulator solves the incompressible Navier-Stokes equations:
          </p>
          <div className="mb-3">
            <Math block>
              {
                "\\frac{\\partial \\mathbf{u}}{\\partial t} + (\\mathbf{u} \\cdot \\nabla)\\mathbf{u} = -\\frac{\\nabla p}{\\rho} + \\nu \\nabla^2 \\mathbf{u} + \\mathbf{f}"
              }
            </Math>
            <Math block>{"\\nabla \\cdot \\mathbf{u} = 0"}</Math>
          </div>
          <p className="mb-2">Where:</p>
          <ul className="list-disc pl-5 space-y-1 text-sm">
            <li>
              <Math>{"\\mathbf{u}"}</Math> = velocity field (u_x, u_y)
            </li>
            <li>
              <Math>{"p"}</Math> = pressure field
            </li>
            <li>
              <Math>{"\\rho"}</Math> = fluid density
            </li>
            <li>
              <Math>{"\\nu"}</Math> = kinematic viscosity
            </li>
            <li>
              <Math>{"\\mathbf{f}"}</Math> = external forces
            </li>
          </ul>
        </div>
      ),
    },

    // Mathematical Foundation
    {
      category: "mathematics",
      question: "How does CFD transform physics equations into computer code?",
      answer: (
        <div>
          <p className="mb-2">
            Computational Fluid Dynamics transforms continuous partial
            differential equations into discrete algebraic equations that
            computers can solve. The journey involves:
          </p>
          <ol className="list-decimal pl-5 space-y-1 mb-2">
            <li>
              <strong>Physical Laws</strong> → Mathematical equations
              (Navier-Stokes)
            </li>
            <li>
              <strong>Continuous Domain</strong> → Discrete grid (finite
              differences)
            </li>
            <li>
              <strong>Time Derivatives</strong> → Time stepping schemes
            </li>
            <li>
              <strong>Spatial Derivatives</strong> → Finite difference
              approximations
            </li>
            <li>
              <strong>Coupled System</strong> → Operator splitting (solve parts
              separately)
            </li>
          </ol>
        </div>
      ),
    },
    {
      category: "mathematics",
      question: "What is the physical meaning of the Navier-Stokes equations?",
      answer: (
        <div>
          <p className="mb-2">
            <strong>Conservation of Mass (Continuity Equation):</strong>
          </p>
          <div className="mb-2">
            <Math block>
              {
                "\\nabla \\cdot \\mathbf{u} = \\frac{\\partial u}{\\partial x} + \\frac{\\partial v}{\\partial y} = 0"
              }
            </Math>
          </div>
          <p className="mb-3 text-sm">
            <strong>Physical meaning:</strong> "What flows in must flow out" -
            mass is conserved.
          </p>

          <p className="mb-2">
            <strong>Conservation of Momentum (Navier-Stokes):</strong>
          </p>
          <div className="mb-2">
            <Math block>
              {
                "\\frac{\\partial \\mathbf{u}}{\\partial t} + (\\mathbf{u} \\cdot \\nabla)\\mathbf{u} = -\\frac{\\nabla p}{\\rho} + \\nu \\nabla^2 \\mathbf{u} + \\mathbf{f}"
              }
            </Math>
          </div>
          <p className="mb-2">Breaking this down term by term:</p>
          <ul className="list-disc pl-5 space-y-1 text-sm">
            <li>
              <strong>
                <Math>{"\\frac{\\partial \\mathbf{u}}{\\partial t}"}</Math>:
              </strong>{" "}
              Time rate of change of velocity (acceleration)
            </li>
            <li>
              <strong>
                <Math>{"(\\mathbf{u} \\cdot \\nabla)\\mathbf{u}"}</Math>:
              </strong>{" "}
              Convective acceleration (fluid accelerating itself)
            </li>
            <li>
              <strong>
                <Math>{"-\\frac{\\nabla p}{\\rho}"}</Math>:
              </strong>{" "}
              Pressure gradient force per unit mass
            </li>
            <li>
              <strong>ν∇²u:</strong> Viscous diffusion (smoothing effect)
            </li>
            <li>
              <strong>f:</strong> External forces (gravity, etc.)
            </li>
          </ul>
          <p className="mt-2 text-sm">
            <strong>Physical meaning:</strong> Newton's second law for fluid
            parcels.
          </p>
        </div>
      ),
    },
    {
      category: "mathematics",
      question: "How does flow around a cylinder work as an example?",
      answer: (
        <div>
          <p className="mb-2">
            Let's trace how the equations describe flow around a circular
            obstacle:
          </p>
          <div className="space-y-3">
            <div>
              <strong>Step 1: Initial Conditions</strong>
              <div className="bg-gray-700 p-2 rounded text-sm font-mono mt-1">
                t = 0: Uniform flow approaches cylinder
                <br />
                u(x,y,0) = [U₀, 0] for x &lt; 0 (upstream)
                <br />
                u(x,y,0) = [0, 0] at obstacle boundary
              </div>
            </div>
            <div>
              <strong>Step 2: Pressure Buildup</strong>
              <p className="text-sm">
                As fluid approaches the cylinder, it must slow down (no-slip
                condition). Velocity decreases → Pressure increases (Bernoulli's
                principle). High pressure forms at stagnation point (front of
                cylinder).
              </p>
            </div>
            <div>
              <strong>Step 3: Flow Separation</strong>
              <p className="text-sm">
                Fluid flows around the cylinder: Top/bottom velocity increases,
                pressure drops. Rear: adverse pressure gradient may cause flow
                separation. Wake formation: low pressure region behind cylinder.
              </p>
            </div>
          </div>
        </div>
      ),
    },
    {
      category: "mathematics",
      question: "What is operator splitting and why is it used?",
      answer: (
        <div>
          <p className="mb-2">
            Instead of solving the full coupled system simultaneously, we split
            it into simpler sub-problems:
          </p>
          <div className="mb-3">
            <strong>Traditional Approach (Difficult):</strong>
            <div className="bg-gray-700 p-2 rounded text-sm font-mono mt-1 mb-2">
              ∂u/∂t + (u·∇)u + ∇p/ρ = ν∇²u
              <br />
              ∇·u = 0
            </div>
            <p className="text-sm">
              This requires solving a large coupled system - computationally
              expensive!
            </p>
          </div>
          <div>
            <strong>Projection Method (Our Approach):</strong>
            <p className="text-sm mb-2">Split into sequential steps:</p>
            <ol className="list-decimal pl-5 space-y-1 text-sm">
              <li>
                <strong>Advection:</strong> ∂u/∂t + (u·∇)u = 0
              </li>
              <li>
                <strong>Diffusion:</strong> ∂u/∂t = ν∇²u
              </li>
              <li>
                <strong>Projection:</strong> Make velocity field divergence-free
              </li>
            </ol>
            <p className="text-sm mt-2">
              <strong>Why this works:</strong> Each sub-problem is easier to
              solve, and the splitting error is small for small time steps.
            </p>
          </div>
        </div>
      ),
    },

    // Architecture & Implementation
    {
      category: "implementation",
      question: "What is the code architecture and structure?",
      answer: (
        <div>
          <div className="mb-3">
            <strong>Core Components:</strong>
            <div className="bg-gray-700 p-3 rounded text-sm font-mono mt-2">
              app/
              <br />
              ├── page.tsx # Main CFD simulator component
              <br />
              ├── components/
              <br />
              │ ├── EducationalContent.tsx # Educational UI
              <br />
              │ ├── AirfoilAnalyzer.tsx # Airfoil analysis
              <br />
              │ ├── PDEStepVisualizer.tsx # PDE visualization
              <br />
              │ └── SimplifiedPDEVisualizer.tsx # Educational viz
              <br />
              └── globals.css # Styling
            </div>
          </div>
          <div>
            <strong>Data Structures:</strong>
            <div className="bg-gray-700 p-3 rounded text-sm font-mono mt-2">
              interface SimulationState &#123;
              <br />
              &nbsp;&nbsp;velocityX: Float32Array; // X-velocity field
              <br />
              &nbsp;&nbsp;velocityY: Float32Array; // Y-velocity field
              <br />
              &nbsp;&nbsp;pressure: Float32Array; // Pressure field
              <br />
              &nbsp;&nbsp;obstacles: boolean[]; // Obstacle mask
              <br />
              &nbsp;&nbsp;particles: Particle[]; // Lagrangian particles
              <br />
              &#125;
            </div>
          </div>
        </div>
      ),
    },
    {
      category: "implementation",
      question: "How is the computational grid set up?",
      answer: (
        <div>
          <p className="mb-2">
            The simulation uses a <strong>staggered grid</strong> (MAC grid)
            where:
          </p>
          <ul className="list-disc pl-5 space-y-1 mb-2">
            <li>Velocities are stored at cell faces</li>
            <li>Pressure is stored at cell centers</li>
            <li>Grid spacing: dx = dy = CELL_SIZE</li>
          </ul>
          <div className="mb-2">
            <strong>Why Staggered Grid?</strong>
            <div className="bg-gray-700 p-3 rounded text-sm font-mono mt-2">
              Regular Grid Problem:
              <br />
              p₁ ---- u ---- p₂ (pressure and velocity at same points)
              <br />
              Could lead to checkerboard pressure patterns!
              <br />
              <br />
              Staggered Grid Solution:
              <br />
              &nbsp;&nbsp;&nbsp;&nbsp;u_&#123;i+1/2&#125;
              <br />
              p_i ——————— p_&#123;i+1&#125; (pressure at centers, velocity at
              faces)
              <br />
              &nbsp;&nbsp;&nbsp;&nbsp;v_&#123;i,j+1/2&#125;
            </div>
          </div>
        </div>
      ),
    },

    // Numerical Methods
    {
      category: "numerical",
      question:
        "How does viscous diffusion work and what is its physical meaning?",
      answer: (
        <div>
          <p className="mb-2">
            <strong>Mathematical Problem:</strong>
          </p>
          <div className="mb-2">
            <Math block>
              {
                "\\frac{\\partial \\mathbf{u}}{\\partial t} = \\nu \\nabla^2 \\mathbf{u}"
              }
            </Math>
          </div>
          <p className="mb-2">
            This is the <strong>diffusion equation</strong> - identical to heat
            conduction!
          </p>

          <p className="mb-2">
            <strong>Physical Interpretation:</strong>
          </p>
          <ul className="list-disc pl-5 space-y-1 mb-2 text-sm">
            <li>
              <strong>High viscosity:</strong> Fluid acts like honey -
              velocities smooth out quickly
            </li>
            <li>
              <strong>Low viscosity:</strong> Fluid acts like water - sharp
              velocity gradients persist
            </li>
            <li>
              <strong>Boundary layers:</strong> Near walls, viscosity creates
              smooth velocity transition
            </li>
          </ul>

          <p className="mb-2">
            <strong>Numerical Approach:</strong>
          </p>
          <div className="mb-2">
            <Math block>
              {
                "\\nabla^2 u \\approx \\frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j}}{\\Delta x^2}"
              }
            </Math>
          </div>
          <p className="text-sm">
            <strong>Stability Condition:</strong>{" "}
            <Math>
              {"\\alpha = \\frac{\\nu \\Delta t}{\\Delta x^2} \\leq 0.25"}
            </Math>{" "}
            (for 2D explicit scheme)
          </p>
        </div>
      ),
    },
    {
      category: "numerical",
      question: "How does pressure projection enforce incompressibility?",
      answer: (
        <div>
          <p className="mb-2">
            <strong>Physical Law:</strong> Mass conservation for incompressible
            flow
          </p>
          <div className="mb-2">
            <Math block>
              {
                "\\nabla \\cdot \\mathbf{u} = \\frac{\\partial u}{\\partial x} + \\frac{\\partial v}{\\partial y} = 0"
              }
            </Math>
          </div>

          <p className="mb-2">
            <strong>Physical meaning:</strong>
          </p>
          <ul className="list-disc pl-5 space-y-1 mb-2 text-sm">
            <li>Divergence &gt; 0: Fluid expanding (unphysical for liquids)</li>
            <li>
              Divergence &lt; 0: Fluid compressing (unphysical for liquids)
            </li>
            <li>Divergence = 0: Mass conserved (physical requirement)</li>
          </ul>

          <p className="mb-2">
            <strong>Helmholtz-Hodge Decomposition:</strong>
          </p>
          <div className="mb-2">
            <Math block>
              {"\\mathbf{u} = \\mathbf{u}_{solenoidal} + \\nabla \\phi"}
            </Math>
            <p className="text-sm mt-2">where:</p>
            <ul className="list-disc pl-5 space-y-1 text-sm">
              <li>
                <Math>{"\\mathbf{u}_{solenoidal}"}</Math>: divergence-free
                (physical velocity)
              </li>
              <li>
                <Math>{"\\nabla \\phi"}</Math>: curl-free (pressure gradient)
              </li>
            </ul>
          </div>

          <p className="mb-2">
            <strong>Solution Process:</strong>
          </p>
          <ol className="list-decimal pl-5 space-y-1 text-sm">
            <li>
              Solve Poisson equation:{" "}
              <Math>{"\\nabla^2 \\phi = \\nabla \\cdot \\mathbf{u}"}</Math>
            </li>
            <li>
              Subtract gradient:{" "}
              <Math>{"\\mathbf{u}_{new} = \\mathbf{u} - \\nabla \\phi"}</Math>
            </li>
          </ol>
        </div>
      ),
    },
    {
      category: "numerical",
      question: "What is Semi-Lagrangian advection and why is it used?",
      answer: (
        <div>
          <p className="mb-2">
            <strong>Mathematical Problem:</strong>
          </p>
          <div className="bg-gray-700 p-2 rounded text-sm font-mono mb-2">
            ∂u/∂t + (u·∇)u = 0
          </div>
          <p className="mb-2">
            This represents <strong>convection</strong> - how the velocity field
            transports itself.
          </p>

          <div className="mb-3">
            <strong>Traditional Eulerian Problem:</strong>
            <div className="bg-gray-700 p-2 rounded text-sm font-mono mt-1 mb-1">
              ∂u/∂t = -u·∇u
              <br />
              Problem: Can be unstable for large velocities (CFL condition)
            </div>
          </div>

          <div className="mb-3">
            <strong>Semi-Lagrangian Solution:</strong>
            <div className="bg-gray-700 p-2 rounded text-sm font-mono mt-1 mb-1">
              Track particles: dx/dt = u(x,t)
              <br />
              u^&#123;n+1&#125;(x) = u^n(x - u·Δt)
              <br />
              Advantage: Unconditionally stable!
            </div>
          </div>

          <p className="mb-2">
            <strong>Algorithm:</strong>
          </p>
          <ol className="list-decimal pl-5 space-y-1 text-sm">
            <li>
              <strong>Trace backward:</strong> Where was this fluid parcel at
              time n?
            </li>
            <li>
              <strong>Interpolate:</strong> What was the velocity at the
              departure point?
            </li>
          </ol>
        </div>
      ),
    },
    {
      category: "numerical",
      question: "How does bilinear interpolation work?",
      answer: (
        <div>
          <p className="mb-2">
            When particles are traced backward in semi-Lagrangian advection,
            they rarely land exactly on grid points. We need to interpolate
            field values at arbitrary positions.
          </p>

          <p className="mb-2">
            <strong>Method:</strong> Given four surrounding grid points,
            interpolate smoothly:
          </p>
          <div className="bg-gray-700 p-3 rounded text-sm font-mono mb-2">
            &nbsp;&nbsp;&nbsp;&nbsp;v01 ——— v11
            <br />
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|
            <br />
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;P&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;P
            = query point (x,y)
            <br />
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|
            <br />
            &nbsp;&nbsp;&nbsp;&nbsp;v00 ——— v10
          </div>

          <p className="mb-2">
            <strong>Algorithm:</strong>
          </p>
          <ol className="list-decimal pl-5 space-y-1 text-sm">
            <li>Find surrounding grid cell: (x₀, y₀) to (x₁, y₁)</li>
            <li>Compute fractional distances: fx = x - x₀, fy = y - y₀</li>
            <li>
              Bilinear interpolation: v(x,y) = v₀₀(1-fx)(1-fy) + v₁₀·fx(1-fy) +
              v₀₁(1-fx)fy + v₁₁·fx·fy
            </li>
          </ol>
        </div>
      ),
    },

    // Features & Usage
    {
      category: "features",
      question: "What visualization modes are available?",
      answer: (
        <div>
          <div className="space-y-3">
            <div>
              <strong>1. Standard Mode</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Velocity field vectors</li>
                <li>Particle traces with speed/pressure coloring</li>
                <li>Streamlines with flow direction indicators</li>
                <li>Boundary layer thickness visualization</li>
              </ul>
            </div>
            <div>
              <strong>2. Pressure Mode</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Enhanced pressure field with contours</li>
                <li>Pressure gradient indicators</li>
                <li>Minimal velocity vectors for reference</li>
              </ul>
            </div>
            <div>
              <strong>3. Smoke Mode</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Advected scalar field visualization</li>
                <li>Multiple color schemes (thermal, rainbow, plasma)</li>
                <li>Density-based opacity</li>
              </ul>
            </div>
          </div>
        </div>
      ),
    },
    {
      category: "features",
      question: "What interactive features are available?",
      answer: (
        <div>
          <ul className="space-y-2">
            <li>
              <strong>Drawing Tools:</strong> Sketch obstacles with brush or
              Bezier curves
            </li>
            <li>
              <strong>Real-time Analysis:</strong> Click anywhere for local flow
              analysis
            </li>
            <li>
              <strong>Parameter Controls:</strong> Adjust viscosity, velocity,
              density in real-time
            </li>
            <li>
              <strong>Educational Mode:</strong> Step-by-step PDE visualization
            </li>
            <li>
              <strong>Reynolds Number Calculator:</strong> Real-time Re
              calculation with flow regime description
            </li>
            <li>
              <strong>Airfoil Analyzer:</strong> Import and analyze airfoil
              geometries
            </li>
          </ul>
        </div>
      ),
    },
    {
      category: "features",
      question: "How does boundary layer visualization work?",
      answer: (
        <div>
          <p className="mb-2">
            The enhanced boundary layer visualization shows:
          </p>
          <ul className="list-disc pl-5 space-y-1 mb-3">
            <li>
              <strong>Thickness indicators:</strong> Yellow lines showing δ₉₉
              (99% velocity recovery)
            </li>
            <li>
              <strong>Velocity profiles:</strong> Curved lines showing velocity
              distribution near walls
            </li>
          </ul>

          <p className="mb-2">
            <strong>Algorithm:</strong>
          </p>
          <ol className="list-decimal pl-5 space-y-1 text-sm">
            <li>
              Find boundary layer thickness where velocity reaches 99% of free
              stream
            </li>
            <li>Draw thickness indicators at regular intervals</li>
            <li>
              Enhanced viscosity treatment near walls creates smooth velocity
              transitions
            </li>
          </ol>
        </div>
      ),
    },
    {
      category: "features",
      question: "How does the particle system work?",
      answer: (
        <div>
          <p className="mb-2">
            Lagrangian particles provide flow visualization by following fluid
            motion:
          </p>

          <div className="bg-gray-700 p-3 rounded text-sm font-mono mb-2">
            interface Particle &#123;
            <br />
            &nbsp;&nbsp;x: number; // Position X<br />
            &nbsp;&nbsp;y: number; // Position Y<br />
            &nbsp;&nbsp;vx: number; // Velocity X<br />
            &nbsp;&nbsp;vy: number; // Velocity Y<br />
            &nbsp;&nbsp;life: number; // Lifetime [0,1]
            <br />
            &#125;
          </div>

          <p className="mb-2">
            <strong>Process:</strong>
          </p>
          <ol className="list-decimal pl-5 space-y-1 text-sm">
            <li>Generate new particles at inlet</li>
            <li>Get velocity from flow field at particle position</li>
            <li>Update position using interpolated velocity</li>
            <li>Remove particles when they leave domain or expire</li>
          </ol>
        </div>
      ),
    },

    // Building Your Own CFD
    {
      category: "building",
      question: "What do I need to understand before building a CFD simulator?",
      answer: (
        <div>
          <p className="mb-2">
            Before diving into implementation, establish these key concepts:
          </p>
          <ol className="list-decimal pl-5 space-y-2">
            <li>
              <strong>Physical Laws:</strong> Navier-Stokes equations govern
              fluid motion
            </li>
            <li>
              <strong>Numerical Methods:</strong> Transform continuous equations
              into discrete algorithms
            </li>
            <li>
              <strong>Data Structures:</strong> Organize velocity, pressure, and
              geometry data efficiently
            </li>
            <li>
              <strong>Visualization:</strong> Render flow fields for analysis
              and debugging
            </li>
            <li>
              <strong>Stability:</strong> Ensure numerical schemes remain stable
              over time
            </li>
          </ol>
        </div>
      ),
    },
    {
      category: "building",
      question: "How do I set up the computational domain and data structures?",
      answer: (
        <div>
          <p className="mb-2">
            <strong>Physical Considerations:</strong>
          </p>
          <div className="bg-gray-700 p-3 rounded text-sm mb-2">
            Domain size: Balance computational cost vs. flow development
            <br />
            - Too small: Boundary effects dominate
            <br />
            - Too large: Wasted computational resources
            <br />- Rule of thumb: 10-20 characteristic lengths downstream
          </div>

          <p className="mb-2">
            <strong>Grid Resolution:</strong>
          </p>
          <div className="bg-gray-700 p-3 rounded text-sm mb-2">
            Cell size determines accuracy and stability:
            <br />
            - Finer grid: Better resolution, higher computational cost
            <br />
            - Coarser grid: Faster computation, potential numerical diffusion
            <br />- Boundary layers need ~10-20 cells for proper resolution
          </div>

          <p className="mb-2">
            <strong>Memory Layout Strategy:</strong>
          </p>
          <div className="bg-gray-700 p-3 rounded text-sm">
            Use typed arrays for performance:
            <br />
            - Float32Array: 32-bit floats for velocity/pressure
            <br />
            - Uint8Array: Boolean arrays for obstacles (memory efficient)
            <br />- Array indexing: row-major order for cache efficiency
          </div>
        </div>
      ),
    },
    {
      category: "building",
      question: "How do I implement boundary conditions?",
      answer: (
        <div>
          <p className="mb-2">
            <strong>Physical Meaning:</strong>
          </p>
          <ol className="list-decimal pl-5 space-y-1 mb-3 text-sm">
            <li>
              <strong>Inlet:</strong> Where fluid enters (Dirichlet condition: u
              = u_specified)
            </li>
            <li>
              <strong>Outlet:</strong> Where fluid exits (Neumann condition:
              ∂u/∂n = 0)
            </li>
            <li>
              <strong>Walls:</strong> Solid boundaries (No-slip: u = 0)
            </li>
            <li>
              <strong>Symmetry:</strong> Mirror conditions (∂u_n/∂n = 0, u_t =
              0)
            </li>
          </ol>

          <p className="mb-2">
            <strong>Mathematical Implementation:</strong>
          </p>
          <ul className="list-disc pl-5 space-y-1 text-sm">
            <li>
              <strong>Dirichlet:</strong> Value specified directly
            </li>
            <li>
              <strong>Neumann:</strong> Gradient specified (use finite
              differences)
            </li>
            <li>
              <strong>Mixed:</strong> Combination based on physical requirements
            </li>
          </ul>
        </div>
      ),
    },
    {
      category: "building",
      question: "What is the complete time integration algorithm?",
      answer: (
        <div>
          <p className="mb-2">
            The solver uses <strong>operator splitting</strong> with the
            following steps:
          </p>
          <ol className="list-decimal pl-5 space-y-2 text-sm">
            <li>
              <strong>Boundary Conditions:</strong> Apply inlet/outlet/wall
              conditions
            </li>
            <li>
              <strong>Viscous Step:</strong> Solve ∂u/∂t = ν∇²u
            </li>
            <li>
              <strong>Projection:</strong> Enforce incompressibility ∇·u = 0
            </li>
            <li>
              <strong>Advection:</strong> Solve ∂u/∂t + (u·∇)u = 0
            </li>
          </ol>

          <div className="bg-gray-700 p-3 rounded text-sm font-mono mt-2">
            For each time step n → n+1:
            <br />
            1. u* = applyBoundaryConditions(u^n)
            <br />
            2. u** = u* + Δt·ν·∇²u*
            <br />
            3. Solve ∇²φ = ∇·u**/Δt
            <br />
            &nbsp;&nbsp;&nbsp;u*** = u** - Δt·∇φ
            <br />
            4. u^&#123;n+1&#125; = advect(u***, u***)
            <br />
            5. u^&#123;n+1&#125; = applyBoundaryConditions(u^&#123;n+1&#125;)
          </div>
        </div>
      ),
    },

    // Performance & Troubleshooting
    {
      category: "performance",
      question: "How can I optimize CFD simulation performance?",
      answer: (
        <div>
          <div className="space-y-3">
            <div>
              <strong>1. Memory Management</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Use Float32Array for large numerical arrays</li>
                <li>Pool particle objects to avoid garbage collection</li>
                <li>Limit maximum particle count</li>
              </ul>
            </div>
            <div>
              <strong>2. Computational Efficiency</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Skip calculations in obstacle cells</li>
                <li>Use appropriate iteration counts for pressure solver</li>
                <li>Implement adaptive time stepping for stability</li>
              </ul>
            </div>
            <div>
              <strong>3. Rendering Optimization</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Use requestAnimationFrame for smooth animation</li>
                <li>Skip expensive operations when not visible</li>
                <li>Implement level-of-detail for particle rendering</li>
              </ul>
            </div>
          </div>
        </div>
      ),
    },
    {
      category: "performance",
      question: "What are the stability considerations?",
      answer: (
        <div>
          <div className="bg-gray-700 p-3 rounded text-sm font-mono mb-2">
            // CFL condition for advection
            <br />
            const maxVelocity = Math.max(
            <br />
            &nbsp;&nbsp;...velocityX.map(Math.abs),
            <br />
            &nbsp;&nbsp;...velocityY.map(Math.abs)
            <br />
            );
            <br />
            const maxTimeStep = (0.5 * CELL_SIZE) / maxVelocity;
            <br />
            const dt = Math.min(0.016, maxTimeStep);
            <br />
            <br />
            // Viscous stability constraint
            <br />
            const viscousTimeStep = (0.25 * CELL_SIZE * CELL_SIZE) / viscosity;
            <br />
            const stableDt = Math.min(dt, viscousTimeStep);
          </div>
          <p className="text-sm">
            Time step must satisfy both CFL condition (advection) and viscous
            stability limit (diffusion) to prevent numerical instability.
          </p>
        </div>
      ),
    },
    {
      category: "performance",
      question: "What are common issues and how do I troubleshoot them?",
      answer: (
        <div>
          <div className="space-y-3">
            <div>
              <strong>Common Issues:</strong>
              <ol className="list-decimal pl-5 space-y-1 text-sm mt-1">
                <li>
                  <strong>Simulation Explodes:</strong> Check CFL condition and
                  reduce time step
                </li>
                <li>
                  <strong>Poor Convergence:</strong> Increase pressure iteration
                  count
                </li>
                <li>
                  <strong>Boundary Layer Not Visible:</strong> Increase
                  viscosity or enable BL visualization
                </li>
                <li>
                  <strong>Performance Issues:</strong> Reduce grid resolution or
                  particle count
                </li>
              </ol>
            </div>
            <div>
              <strong>Debugging Tools:</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Add velocity magnitude visualization</li>
                <li>Monitor maximum velocities and pressures</li>
                <li>Implement divergence checking for mass conservation</li>
                <li>Use console logging for critical parameters</li>
              </ul>
            </div>
            <div>
              <strong>Physical Validation:</strong>
              <ul className="list-disc pl-5 space-y-1 text-sm mt-1">
                <li>Verify Reynolds number calculations</li>
                <li>Check mass conservation (∇·u ≈ 0)</li>
                <li>Validate boundary conditions</li>
                <li>
                  Compare results with analytical solutions for simple cases
                </li>
              </ul>
            </div>
          </div>
        </div>
      ),
    },

    // Installation & Setup
    {
      category: "setup",
      question: "How do I install and run this CFD simulator?",
      answer: (
        <div>
          <div className="bg-gray-700 p-3 rounded text-sm font-mono mb-3">
            # Clone repository
            <br />
            git clone &lt;https://github.com/rongeld/cfd-browser&gt;
            <br />
            cd cfd-simulator
            <br />
            <br />
            # Install dependencies
            <br />
            npm install
            <br />
            <br />
            # Run development server
            <br />
            npm run dev
            <br />
            <br />
            # Build for production
            <br />
            npm run build
            <br />
            npm start
          </div>
          <div>
            <strong>Dependencies:</strong>
            <div className="bg-gray-700 p-3 rounded text-sm font-mono mt-2">
              &#123;
              <br />
              &nbsp;&nbsp;"next": "^14.0.0",
              <br />
              &nbsp;&nbsp;"react": "^18.0.0",
              <br />
              &nbsp;&nbsp;"typescript": "^5.0.0",
              <br />
              &nbsp;&nbsp;"lucide-react": "^0.263.1",
              <br />
              &nbsp;&nbsp;"recharts": "^2.8.0"
              <br />
              &#125;
            </div>
          </div>
        </div>
      ),
    },
  ];

  const categories = [
    { id: "all", name: "All Questions" },
    { id: "overview", name: "Overview & Basics" },
    { id: "mathematics", name: "Mathematical Foundation" },
    { id: "implementation", name: "Architecture & Implementation" },
    { id: "numerical", name: "Numerical Methods" },
    { id: "features", name: "Features & Visualization" },
    { id: "building", name: "Building Your Own CFD" },
    { id: "performance", name: "Performance & Troubleshooting" },
    { id: "setup", name: "Installation & Setup" },
  ];

  const filteredFAQ =
    activeCategory === "all"
      ? faqData
      : faqData.filter((item) => item.category === activeCategory);

  return (
    <div className={`bg-gray-900 text-white ${className}`}>
      <div className="max-w-6xl mx-auto px-6 py-8">
        {/* Header */}
        <div className="text-center mb-8">
          <h1 className="text-4xl font-bold mb-4">
            Frequently Asked Questions
          </h1>
          <p className="text-gray-400 text-lg">
            Everything you need to know about CFD simulation and this
            application
          </p>
        </div>

        {/* Category Filter */}
        <div className="mb-8">
          <div className="flex flex-wrap justify-center gap-2">
            {categories.map((category) => (
              <button
                key={category.id}
                onClick={() => setActiveCategory(category.id)}
                className={`px-4 py-2 rounded-lg text-sm transition-colors ${
                  activeCategory === category.id
                    ? "bg-blue-600 text-white"
                    : "bg-gray-700 text-gray-300 hover:bg-gray-600"
                }`}
              >
                {category.name}
              </button>
            ))}
          </div>
        </div>

        {/* FAQ Items */}
        <div className="space-y-4">
          {filteredFAQ.map((item, index) => (
            <div
              key={index}
              className="bg-gray-800 rounded-lg border border-gray-700"
            >
              <button
                onClick={() => toggleItem(index)}
                className="w-full px-6 py-4 text-left flex items-center justify-between hover:bg-gray-750 transition-colors"
              >
                <h3 className="text-lg font-semibold pr-4">{item.question}</h3>
                {openItems.has(index) ? (
                  <ChevronDown className="w-5 h-5 text-blue-400 flex-shrink-0" />
                ) : (
                  <ChevronRight className="w-5 h-5 text-gray-400 flex-shrink-0" />
                )}
              </button>
              {openItems.has(index) && (
                <div className="px-6 pb-4 text-gray-300 border-t border-gray-700">
                  <div className="pt-4">{item.answer}</div>
                </div>
              )}
            </div>
          ))}
        </div>

        {/* Additional Resources */}
        <div className="mt-12 text-center">
          <div className="bg-gray-800 rounded-lg p-6 border border-gray-700">
            <h3 className="text-xl font-semibold mb-3">
              Still have questions?
            </h3>
            <p className="text-gray-400 mb-4">
              For more detailed technical information, check out the
              comprehensive documentation in the "Learn CFD" tab.
            </p>
            <div className="text-sm text-gray-500">
              This simulator is open source and educational. Feel free to
              explore the code and build upon it!
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default FAQSection;
