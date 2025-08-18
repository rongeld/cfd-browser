"use client";
import React, { useState, useRef, useEffect } from "react";
import {
  ChevronRight,
  ChevronLeft,
  Code,
  Play,
  RotateCcw,
  Eye,
  EyeOff,
  BookOpen,
  Zap,
  Wind,
  Calculator,
} from "lucide-react";

interface Step {
  id: number;
  title: string;
  description: string;
  physicsExplanation: string;
  codeChanges: {
    additions: string[];
    explanations: string[];
  };
  executionResult: {
    description: string;
    visualEffect: string;
  };
}

const InteractiveCFDBuilder: React.FC = () => {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const animationRef = useRef<number | null>(null);
  const codeScrollRef = useRef<HTMLDivElement | null>(null);
  const [currentStep, setCurrentStep] = useState(0);
  const [isRunning, setIsRunning] = useState(false);
  const [showCode, setShowCode] = useState(true);
  const [showPhysics, setShowPhysics] = useState(true);
  const [isScrolling, setIsScrolling] = useState(false);

  // Build cumulative code up to current step
  const getCumulativeCode = () => {
    const cumulativeLines: Array<{
      text: string;
      stepId: number;
      isNew: boolean;
      isStepSeparator?: boolean;
    }> = [];

    for (let i = 0; i <= currentStep; i++) {
      const step = steps[i];
      step.codeChanges.additions.forEach((line) => {
        cumulativeLines.push({
          text: line,
          stepId: i,
          isNew: i === currentStep,
          isStepSeparator: false,
        });
      });

      // Add separator between steps (except for last step)
      if (i < currentStep) {
        cumulativeLines.push({
          text: "",
          stepId: i,
          isNew: false,
          isStepSeparator: false,
        });
        cumulativeLines.push({
          text: `// ==================== Step ${i + 2}: ${
            steps[i + 1]?.title
          } ====================`,
          stepId: i,
          isNew: false,
          isStepSeparator: true,
        });
        cumulativeLines.push({
          text: "",
          stepId: i,
          isNew: false,
          isStepSeparator: false,
        });
      }
    }

    return cumulativeLines;
  };

  // Scroll to newly added code when step changes
  const scrollToNewCode = () => {
    if (!codeScrollRef.current || !showCode) return;

    const container = codeScrollRef.current;

    // Wait for DOM to fully update
    setTimeout(() => {
      // Try to find the first highlighted element with actual content
      const highlightedElements = container.querySelectorAll(
        '[data-is-new="true"]'
      );
      let targetElement = null;

      // Find the first highlighted element that has non-empty text content
      for (let i = 0; i < highlightedElements.length; i++) {
        const element = highlightedElements[i] as HTMLElement;
        const textContent = element.textContent?.trim() || "";
        if (textContent !== "" && !textContent.startsWith("//")) {
          targetElement = element;
          break;
        }
      }

      if (targetElement) {
        setIsScrolling(true);

        // Calculate scroll position - we want to show the new content near the top
        const elementOffsetTop = (targetElement as HTMLElement).offsetTop;

        // Target position: show new code about 60px from the top of visible area
        const targetScrollTop = elementOffsetTop - 60;

        container.scrollTo({
          top: Math.max(0, targetScrollTop),
          behavior: "smooth",
        });

        // Hide scrolling indicator after animation
        setTimeout(() => {
          setIsScrolling(false);
        }, 1200);
      } else {
        // Fallback: scroll toward the end where new content appears
        const maxScroll = container.scrollHeight - container.clientHeight;
        if (maxScroll > 0) {
          setIsScrolling(true);
          container.scrollTo({
            top: Math.max(0, maxScroll - 50), // Leave some space at bottom
            behavior: "smooth",
          });
          setTimeout(() => {
            setIsScrolling(false);
          }, 1200);
        }
      }
    }, 200); // Longer delay to ensure rendering is complete
  };

  // Trigger scroll when step changes
  useEffect(() => {
    if (showCode) {
      // Longer delay to ensure the DOM has fully updated with new content
      const timer = setTimeout(() => {
        scrollToNewCode();
      }, 150);

      return () => clearTimeout(timer);
    }
  }, [currentStep, showCode]);

  const steps: Step[] = [
    {
      id: 0,
      title: "Initialize Canvas and Grid",
      description:
        "Set up our computational domain - the space where fluid will flow",
      physicsExplanation:
        "In CFD, we discretize space into a computational grid or mesh. Each cell represents a small volume of fluid where we'll calculate properties like velocity, pressure, and density. Think of it like dividing a river into thousands of tiny boxes - each box tells us what the water is doing there. The finer the grid, the more accurate our simulation, but also the more computationally expensive.",
      codeChanges: {
        additions: [
          "const canvas = document.getElementById('cfd-canvas');",
          "const ctx = canvas.getContext('2d');",
          "const GRID_WIDTH = 120;   // Number of cells horizontally",
          "const GRID_HEIGHT = 60;   // Number of cells vertically",
          "const CELL_SIZE = 6;      // Size of each cell in pixels",
          "",
          "// Initialize fluid properties for each grid cell",
          "let velocityX = new Float32Array(GRID_WIDTH * GRID_HEIGHT);",
          "let velocityY = new Float32Array(GRID_WIDTH * GRID_HEIGHT);",
          "let pressure = new Float32Array(GRID_WIDTH * GRID_HEIGHT);",
          "let density = new Float32Array(GRID_WIDTH * GRID_HEIGHT);",
          "",
          "// Fill with default values",
          "density.fill(1.0); // Standard air density kg/m³",
        ],
        explanations: [
          "Create 2D rendering context for visualization",
          "Define computational grid dimensions",
          "Use Float32Array for better performance",
          "Initialize velocity components (u, v)",
          "Initialize pressure and density fields",
          "Set default fluid properties",
        ],
      },
      executionResult: {
        description:
          "Computational grid established with 7,200 cells ready for simulation",
        visualEffect: "grid",
      },
    },
    {
      id: 1,
      title: "Set Initial Flow Conditions",
      description: "Create uniform flow from left to right - our 'wind tunnel'",
      physicsExplanation:
        "We establish initial conditions by setting a uniform velocity field. This represents air or water flowing at constant speed from left to right. The Reynolds number Re = UL/ν determines flow characteristics: low Re means laminar (smooth) flow, high Re leads to turbulence. For our circle, Re ≈ 1000 creates beautiful vortex shedding.",
      codeChanges: {
        additions: [
          "const INLET_VELOCITY = 8.0;    // m/s - inlet flow speed",
          "const REYNOLDS_NUMBER = 1000;   // Dimensionless flow parameter",
          "const KINEMATIC_VISCOSITY = INLET_VELOCITY * OBSTACLE_RADIUS / REYNOLDS_NUMBER;",
          "",
          "// Initialize uniform flow field",
          "function initializeFlow() {",
          "  for (let y = 0; y < GRID_HEIGHT; y++) {",
          "    for (let x = 0; x < GRID_WIDTH; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      velocityX[index] = INLET_VELOCITY;",
          "      velocityY[index] = 0.0;",
          "      pressure[index] = 0.0;",
          "    }",
          "  }",
          "}",
          "",
          "initializeFlow();",
        ],
        explanations: [
          "Set realistic inlet velocity for wind tunnel",
          "Define Reynolds number for flow regime",
          "Calculate viscosity from Re = UL/ν",
          "Loop through all computational cells",
          "Set uniform horizontal velocity",
          "Initialize with zero vertical velocity and pressure",
        ],
      },
      executionResult: {
        description: "Uniform flow field established - air moving at 8 m/s",
        visualEffect: "arrows",
      },
    },
    {
      id: 2,
      title: "Add Circular Obstacle",
      description:
        "Place a cylinder in the flow - the classic 'flow around cylinder' problem",
      physicsExplanation:
        "The circular cylinder is one of the most studied geometries in fluid mechanics. When flow encounters the cylinder, it creates a stagnation point at the front (zero velocity, maximum pressure), accelerated flow around the sides (high velocity, low pressure per Bernoulli's principle), and complex wake dynamics behind it. This is fundamental to understanding lift and drag forces.",
      codeChanges: {
        additions: [
          "const OBSTACLE_X = GRID_WIDTH * 0.25;  // Position at 25% from inlet",
          "const OBSTACLE_Y = GRID_HEIGHT * 0.5;  // Center vertically",
          "const OBSTACLE_RADIUS = 6;             // Cylinder radius",
          "",
          "// Create solid boundary marker",
          "let isSolid = new Boolean(GRID_WIDTH * GRID_HEIGHT).fill(false);",
          "",
          "function createObstacle() {",
          "  for (let y = 0; y < GRID_HEIGHT; y++) {",
          "    for (let x = 0; x < GRID_WIDTH; x++) {",
          "      const dx = x - OBSTACLE_X;",
          "      const dy = y - OBSTACLE_Y;",
          "      const distance = Math.sqrt(dx * dx + dy * dy);",
          "      ",
          "      if (distance <= OBSTACLE_RADIUS) {",
          "        const index = y * GRID_WIDTH + x;",
          "        isSolid[index] = true;",
          "        // Enforce no-slip boundary condition",
          "        velocityX[index] = 0.0;",
          "        velocityY[index] = 0.0;",
          "      }",
          "    }",
          "  }",
          "}",
          "",
          "createObstacle();",
        ],
        explanations: [
          "Position obstacle for optimal flow development",
          "Create boolean array for solid/fluid identification",
          "Calculate Euclidean distance to cylinder center",
          "Mark interior cells as solid boundaries",
          "Apply no-slip condition: fluid sticks to walls",
          "Zero velocity enforced at all solid boundaries",
        ],
      },
      executionResult: {
        description:
          "Circular cylinder obstacle added - flow will now have to navigate around it",
        visualEffect: "circle",
      },
    },
    {
      id: 3,
      title: "Implement Mass Conservation",
      description:
        "Ensure fluid is incompressible - what goes in must come out",
      physicsExplanation:
        "The continuity equation ∇·u = 0 states that for incompressible flow, the divergence of velocity must be zero. This means mass is conserved - if fluid flows into a region, the same amount must flow out. We enforce this through the pressure projection method, solving ∇²p = ∇·u* to find pressure corrections that make the velocity field divergence-free.",
      codeChanges: {
        additions: [
          "function calculateDivergence(velX, velY) {",
          "  const divergence = new Float32Array(GRID_WIDTH * GRID_HEIGHT);",
          "  ",
          "  for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "    for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      if (isSolid[index]) continue;",
          "      ",
          "      // Central difference approximation of divergence",
          "      const dudx = 0.5 * (velX[index + 1] - velX[index - 1]);",
          "      const dvdy = 0.5 * (velY[index + GRID_WIDTH] - velY[index - GRID_WIDTH]);",
          "      ",
          "      divergence[index] = dudx + dvdy;",
          "    }",
          "  }",
          "  ",
          "  return divergence;",
          "}",
          "",
          "// Check mass conservation quality",
          "function getMassConservationError() {",
          "  const div = calculateDivergence(velocityX, velocityY);",
          "  let maxDivergence = 0;",
          "  for (let i = 0; i < div.length; i++) {",
          "    maxDivergence = Math.max(maxDivergence, Math.abs(div[i]));",
          "  }",
          "  return maxDivergence;",
          "}",
        ],
        explanations: [
          "Compute velocity divergence using finite differences",
          "Central difference scheme for accuracy",
          "Skip solid cells in calculation",
          "Monitor mass conservation quality",
          "Track maximum divergence as error metric",
        ],
      },
      executionResult: {
        description:
          "Mass conservation monitoring implemented - tracking incompressibility",
        visualEffect: "divergence",
      },
    },
    {
      id: 4,
      title: "Add Momentum Advection",
      description:
        "Transport momentum with the flow - make fluid carry its motion",
      physicsExplanation:
        "The advection term u·∇u in Navier-Stokes represents momentum transport. When fast-moving fluid encounters slow fluid, momentum is transferred. This nonlinear term is responsible for many interesting flow phenomena like vortex formation and turbulence. We use upwind differencing to maintain stability in high-speed flows.",
      codeChanges: {
        additions: [
          "function advectVelocity(velField, velX, velY, dt) {",
          "  const newVel = new Float32Array(velField.length);",
          "  ",
          "  for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "    for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      if (isSolid[index]) continue;",
          "      ",
          "      const u = velX[index];",
          "      const v = velY[index];",
          "      ",
          "      // Semi-Lagrangian advection: trace backwards",
          "      const prevX = x - u * dt;",
          "      const prevY = y - v * dt;",
          "      ",
          "      // Bilinear interpolation for smooth transport",
          "      newVel[index] = interpolateField(velField, prevX, prevY);",
          "    }",
          "  }",
          "  ",
          "  return newVel;",
          "}",
          "",
          "function interpolateField(field, x, y) {",
          "  const x0 = Math.floor(x);",
          "  const y0 = Math.floor(y);",
          "  const fx = x - x0;",
          "  const fy = y - y0;",
          "  ",
          "  if (x0 < 0 || x0 >= GRID_WIDTH-1 || y0 < 0 || y0 >= GRID_HEIGHT-1) {",
          "    return 0;",
          "  }",
          "  ",
          "  const i00 = y0 * GRID_WIDTH + x0;",
          "  const i10 = y0 * GRID_WIDTH + (x0 + 1);",
          "  const i01 = (y0 + 1) * GRID_WIDTH + x0;",
          "  const i11 = (y0 + 1) * GRID_WIDTH + (x0 + 1);",
          "  ",
          "  return field[i00] * (1-fx) * (1-fy) +",
          "         field[i10] * fx * (1-fy) +",
          "         field[i01] * (1-fx) * fy +",
          "         field[i11] * fx * fy;",
          "}",
        ],
        explanations: [
          "Semi-Lagrangian method for stable advection",
          "Trace particle paths backwards in time",
          "Bilinear interpolation preserves smoothness",
          "Handle boundary conditions properly",
          "Maintain momentum conservation during transport",
        ],
      },
      executionResult: {
        description:
          "Momentum advection active - flow patterns now develop naturally",
        visualEffect: "advection",
      },
    },
    {
      id: 3,
      title: "Implement Advection",
      description:
        "Make fluid particles move with the flow - transport properties downstream",
      physicsExplanation:
        "Advection describes how fluid properties are transported by the flow. If you put dye in water, it moves downstream - that's advection. Mathematically, it's ∂φ/∂t + u·∇φ = 0, where φ is any transported quantity and u is velocity.",
      codeChanges: {
        additions: [
          "function advect(field, velocityX, velocityY, dt) {",
          "  const newField = [...field];",
          "  ",
          "  for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "    for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      ",
          "      if (obstacle[index]) continue;",
          "      ",
          "      // Trace particle backwards in time",
          "      const prevX = x - velocityX[index] * dt;",
          "      const prevY = y - velocityY[index] * dt;",
          "      ",
          "      // Bilinear interpolation",
          "      const x0 = Math.floor(prevX);",
          "      const y0 = Math.floor(prevY);",
          "      const fx = prevX - x0;",
          "      const fy = prevY - y0;",
          "      ",
          "      if (x0 >= 0 && x0 < GRID_WIDTH-1 && y0 >= 0 && y0 < GRID_HEIGHT-1) {",
          "        const i00 = y0 * GRID_WIDTH + x0;",
          "        const i10 = y0 * GRID_WIDTH + (x0 + 1);",
          "        const i01 = (y0 + 1) * GRID_WIDTH + x0;",
          "        const i11 = (y0 + 1) * GRID_WIDTH + (x0 + 1);",
          "        ",
          "        newField[index] = ",
          "          field[i00] * (1 - fx) * (1 - fy) +",
          "          field[i10] * fx * (1 - fy) +",
          "          field[i01] * (1 - fx) * fy +",
          "          field[i11] * fx * fy;",
          "      }",
          "    }",
          "  }",
          "  ",
          "  return newField;",
          "}",
        ],
        explanations: [
          "Semi-Lagrangian advection scheme",
          "Trace particle position backwards in time",
          "Use bilinear interpolation for smooth transport",
          "Maintain conservation properties",
        ],
      },
      executionResult: {
        description: "Fluid properties now transport with the flow",
        visualEffect: "advection",
      },
    },
    {
      id: 4,
      title: "Add Pressure Projection",
      description:
        "Enforce incompressibility - make sure fluid doesn't compress or expand",
      physicsExplanation:
        "Incompressible flow means ∇·u = 0 (divergence-free velocity field). We solve a pressure Poisson equation to find pressure corrections that make the velocity field divergence-free. This is the heart of the projection method used in most CFD codes.",
      codeChanges: {
        additions: [
          "function project(velocityX, velocityY) {",
          "  const divergence = new Array(GRID_WIDTH * GRID_HEIGHT).fill(0);",
          "  const pressureCorrection = new Array(GRID_WIDTH * GRID_HEIGHT).fill(0);",
          "  ",
          "  // Calculate divergence of velocity field",
          "  for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "    for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      if (obstacle[index]) continue;",
          "      ",
          "      divergence[index] = 0.5 * (",
          "        velocityX[index + 1] - velocityX[index - 1] +",
          "        velocityY[index + GRID_WIDTH] - velocityY[index - GRID_WIDTH]",
          "      );",
          "    }",
          "  }",
          "  ",
          "  // Solve pressure Poisson equation (Jacobi iterations)",
          "  for (let iter = 0; iter < 20; iter++) {",
          "    const newPressure = [...pressureCorrection];",
          "    for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "      for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "        const index = y * GRID_WIDTH + x;",
          "        if (obstacle[index]) continue;",
          "        ",
          "        newPressure[index] = 0.25 * (",
          "          pressureCorrection[index - 1] + pressureCorrection[index + 1] +",
          "          pressureCorrection[index - GRID_WIDTH] + pressureCorrection[index + GRID_WIDTH] -",
          "          divergence[index]",
          "        );",
          "      }",
          "    }",
          "    pressureCorrection = newPressure;",
          "  }",
          "  ",
          "  // Apply pressure correction to velocity",
          "  for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "    for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      if (obstacle[index]) continue;",
          "      ",
          "      velocityX[index] -= 0.5 * (pressureCorrection[index + 1] - pressureCorrection[index - 1]);",
          "      velocityY[index] -= 0.5 * (pressureCorrection[index + GRID_WIDTH] - pressureCorrection[index - GRID_WIDTH]);",
          "    }",
          "  }",
          "}",
        ],
        explanations: [
          "Calculate velocity divergence (compressibility)",
          "Solve Poisson equation for pressure",
          "Apply pressure gradient to make flow incompressible",
          "Ensures mass conservation",
        ],
      },
      executionResult: {
        description: "Flow is now incompressible - mass is conserved",
        visualEffect: "projection",
      },
    },
    {
      id: 5,
      title: "Apply Boundary Conditions",
      description: "Set proper conditions at walls and inlets/outlets",
      physicsExplanation:
        "Boundary conditions define how fluid behaves at domain edges. No-slip means fluid sticks to solid walls (velocity = 0). At the inlet we specify velocity, at the outlet we let flow exit freely. These conditions are crucial for realistic physics.",
      codeChanges: {
        additions: [
          "function applyBoundaryConditions() {",
          "  // Inlet boundary (left side) - specify velocity",
          "  for (let y = 0; y < GRID_HEIGHT; y++) {",
          "    const index = y * GRID_WIDTH + 0;",
          "    velocityX[index] = INITIAL_VELOCITY;",
          "    velocityY[index] = 0;",
          "  }",
          "  ",
          "  // Outlet boundary (right side) - free outflow",
          "  for (let y = 0; y < GRID_HEIGHT; y++) {",
          "    const index = y * GRID_WIDTH + (GRID_WIDTH - 1);",
          "    velocityX[index] = velocityX[index - 1];",
          "    velocityY[index] = velocityY[index - 1];",
          "  }",
          "  ",
          "  // Top and bottom walls - no-slip condition",
          "  for (let x = 0; x < GRID_WIDTH; x++) {",
          "    // Top wall",
          "    velocityX[x] = 0;",
          "    velocityY[x] = 0;",
          "    // Bottom wall",
          "    const bottomIndex = (GRID_HEIGHT - 1) * GRID_WIDTH + x;",
          "    velocityX[bottomIndex] = 0;",
          "    velocityY[bottomIndex] = 0;",
          "  }",
          "  ",
          "  // Obstacle boundary - no-slip",
          "  for (let y = 0; y < GRID_HEIGHT; y++) {",
          "    for (let x = 0; x < GRID_WIDTH; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      if (obstacle[index]) {",
          "        velocityX[index] = 0;",
          "        velocityY[index] = 0;",
          "      }",
          "    }",
          "  }",
          "}",
        ],
        explanations: [
          "Inlet: constant velocity inflow",
          "Outlet: zero gradient (natural outflow)",
          "Walls: no-slip condition (velocity = 0)",
          "Obstacle: no-slip at solid boundaries",
        ],
      },
      executionResult: {
        description: "Proper boundary conditions applied",
        visualEffect: "boundaries",
      },
    },
    {
      id: 6,
      title: "Add Viscous Diffusion",
      description: "Include viscosity effects - make fluid 'sticky'",
      physicsExplanation:
        "Viscosity causes momentum to diffuse through the fluid, smoothing out sharp velocity gradients. This is the ν∇²u term in the Navier-Stokes equations. Higher viscosity means more diffusion and smoother flow patterns.",
      codeChanges: {
        additions: [
          "const VISCOSITY = 0.1; // Kinematic viscosity",
          "",
          "function diffuse(field, viscosity, dt) {",
          "  const alpha = dt * viscosity * GRID_WIDTH * GRID_HEIGHT;",
          "  const beta = 1 + 4 * alpha;",
          "  ",
          "  // Gauss-Seidel iterations for diffusion equation",
          "  for (let iter = 0; iter < 20; iter++) {",
          "    for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "      for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "        const index = y * GRID_WIDTH + x;",
          "        if (obstacle[index]) continue;",
          "        ",
          "        field[index] = (",
          "          field[index] + alpha * (",
          "            field[index - 1] + field[index + 1] +",
          "            field[index - GRID_WIDTH] + field[index + GRID_WIDTH]",
          "          )",
          "        )",
          "      }",
          "    }",
          "  }",
          "}",
        ],
        explanations: [
          "Implicit diffusion for stability",
          "Gauss-Seidel solver for efficiency",
          "Smooths velocity gradients",
          "Represents molecular viscosity effects",
        ],
      },
      executionResult: {
        description: "Viscous effects added - smoother flow",
        visualEffect: "diffusion",
      },
    },
    {
      id: 7,
      title: "Complete Time Integration",
      description: "Put it all together - the main simulation loop",
      physicsExplanation:
        "We use operator splitting: advect velocities, add forces, diffuse, then project to maintain incompressibility. This is the fractional step method, widely used in CFD. Each step advances the solution forward in time by dt.",
      codeChanges: {
        additions: [
          "const dt = 0.1; // Time step",
          "",
          "function simulationStep() {",
          "  // Step 1: Advection",
          "  velocityX = advect(velocityX, velocityX, velocityY, dt);",
          "  velocityY = advect(velocityY, velocityX, velocityY, dt);",
          "  ",
          "  // Step 2: Apply boundary conditions",
          "  applyBoundaryConditions();",
          "  ",
          "  // Step 3: Viscous diffusion",
          "  diffuse(velocityX, VISCOSITY, dt);",
          "  diffuse(velocityY, VISCOSITY, dt);",
          "  ",
          "  // Step 4: Apply boundary conditions again",
          "  applyBoundaryConditions();",
          "  ",
          "  // Step 5: Pressure projection (enforce incompressibility)",
          "  project(velocityX, velocityY);",
          "  ",
          "  // Step 6: Final boundary conditions",
          "  applyBoundaryConditions();",
          "}",
          "",
          "function animate() {",
          "  simulationStep();",
          "  render();",
          "  requestAnimationFrame(animate);",
          "}",
          "",
          "// Start simulation",
          "animate();",
        ],
        explanations: [
          "Fractional step method",
          "Operator splitting for stability",
          "Each step preserves different physics",
          "Animation loop for real-time visualization",
        ],
      },
      executionResult: {
        description: "Complete CFD simulation running!",
        visualEffect: "simulation",
      },
    },
    {
      id: 8,
      title: "Add Vorticity Visualization",
      description:
        "See the spinning motion in the fluid - where turbulence begins",
      physicsExplanation:
        "Vorticity ω = ∇×u measures rotation in the flow field. It's a fundamental quantity in fluid mechanics - high vorticity regions indicate strong shear and potential instability. Around cylinders, vortices shed periodically (Kármán vortex street), creating alternating lift forces. This is why flags flutter and why bridges can oscillate in wind.",
      codeChanges: {
        additions: [
          "function calculateVorticity(velocityX, velocityY) {",
          "  const vorticity = new Float32Array(GRID_WIDTH * GRID_HEIGHT);",
          "  ",
          "  for (let y = 1; y < GRID_HEIGHT - 1; y++) {",
          "    for (let x = 1; x < GRID_WIDTH - 1; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      if (isSolid[index]) continue;",
          "      ",
          "      // Vorticity: ω = ∂v/∂x - ∂u/∂y",
          "      const dvdx = 0.5 * (velocityY[index + 1] - velocityY[index - 1]);",
          "      const dudy = 0.5 * (velocityX[index + GRID_WIDTH] - velocityX[index - GRID_WIDTH]);",
          "      vorticity[index] = dvdx - dudy;",
          "    }",
          "  }",
          "  ",
          "  return vorticity;",
          "}",
          "",
          "// Enhanced visualization with vorticity contours",
          "function renderVorticityField() {",
          "  const vorticity = calculateVorticity(velocityX, velocityY);",
          "  ",
          "  for (let y = 0; y < GRID_HEIGHT; y++) {",
          "    for (let x = 0; x < GRID_WIDTH; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      if (isSolid[index]) continue;",
          "      ",
          "      const omega = vorticity[index];",
          "      const intensity = Math.min(1, Math.abs(omega) * 10);",
          "      ",
          "      // Color code: red for clockwise, blue for counter-clockwise",
          "      if (omega > 0.01) {",
          "        ctx.fillStyle = `rgba(59, 130, 246, ${intensity * 0.6})`;",
          "      } else if (omega < -0.01) {",
          "        ctx.fillStyle = `rgba(239, 68, 68, ${intensity * 0.6})`;",
          "      } else {",
          "        continue;",
          "      }",
          "      ",
          "      const cellX = x * CELL_SIZE;",
          "      const cellY = y * CELL_SIZE;",
          "      ctx.fillRect(cellX, cellY, CELL_SIZE, CELL_SIZE);",
          "    }",
          "  }",
          "}",
        ],
        explanations: [
          "Calculate curl of velocity field",
          "Identify regions of fluid rotation",
          "Color-code vorticity magnitude and sign",
          "Visualize vortex shedding behind cylinder",
          "Foundation for turbulence analysis",
        ],
      },
      executionResult: {
        description:
          "Vorticity field visible - see the fluid rotation patterns!",
        visualEffect: "vorticity",
      },
    },
    {
      id: 9,
      title: "Engineering Analysis Tools",
      description: "Calculate drag, lift, and performance metrics",
      physicsExplanation:
        "Engineering CFD is about extracting useful data. Drag coefficient Cd relates force to dynamic pressure: Drag = ½ρU²ACd. For cylinders, Cd ≈ 1.2 at Re~1000. The Strouhal number St = fD/U characterizes vortex shedding frequency - critical for avoiding resonance in structures. These dimensionless numbers allow scaling from lab models to full-size applications.",
      codeChanges: {
        additions: [
          "function calculateForces() {",
          "  let dragForce = 0;",
          "  let liftForce = 0;",
          "  ",
          "  // Integrate pressure and viscous forces around cylinder",
          "  for (let y = 0; y < GRID_HEIGHT; y++) {",
          "    for (let x = 0; x < GRID_WIDTH; x++) {",
          "      const index = y * GRID_WIDTH + x;",
          "      ",
          "      if (isObstacleBoundary(x, y)) {",
          "        const dx = x - OBSTACLE_X;",
          "        const dy = y - OBSTACLE_Y;",
          "        const distance = Math.sqrt(dx * dx + dy * dy);",
          "        ",
          "        // Unit normal vector pointing outward",
          "        const nx = dx / distance;",
          "        const ny = dy / distance;",
          "        ",
          "        // Estimate pressure from Bernoulli equation",
          "        const neighborIndex = findNearestFluidCell(x, y);",
          "        const velMag = Math.sqrt(",
          "          velocityX[neighborIndex] ** 2 + velocityY[neighborIndex] ** 2",
          "        );",
          "        const dynamicPressure = 0.5 * velMag * velMag;",
          "        const pressureCoeff = 1 - dynamicPressure / (0.5 * INLET_VELOCITY ** 2);",
          "        ",
          "        // Force = pressure × area",
          "        const pressureForce = pressureCoeff * CELL_SIZE;",
          "        dragForce += pressureForce * nx;",
          "        liftForce += pressureForce * ny;",
          "      }",
          "    }",
          "  }",
          "  ",
          "  // Non-dimensionalize to get coefficients",
          "  const dynamicPressureRef = 0.5 * INLET_VELOCITY ** 2;",
          "  const referenceArea = 2 * OBSTACLE_RADIUS * CELL_SIZE;",
          "  ",
          "  const dragCoefficient = dragForce / (dynamicPressureRef * referenceArea);",
          "  const liftCoefficient = liftForce / (dynamicPressureRef * referenceArea);",
          "  ",
          "  return { dragCoeff: dragCoefficient, liftCoeff: liftCoefficient };",
          "}",
          "",
          "function calculateReynoldsNumber() {",
          "  const characteristicLength = 2 * OBSTACLE_RADIUS * CELL_SIZE;",
          "  const kinematicViscosity = 0.0001; // m²/s for our simulation",
          "  return (INLET_VELOCITY * characteristicLength) / kinematicViscosity;",
          "}",
          "",
          "function isObstacleBoundary(x, y) {",
          "  const dx = x - OBSTACLE_X;",
          "  const dy = y - OBSTACLE_Y;",
          "  const distance = Math.sqrt(dx * dx + dy * dy);",
          "  return Math.abs(distance - OBSTACLE_RADIUS) < 1.0;",
          "}",
          "",
          "function findNearestFluidCell(x, y) {",
          "  // Find nearest non-solid cell for pressure estimation",
          "  for (let r = 1; r < 5; r++) {",
          "    for (let dy = -r; dy <= r; dy++) {",
          "      for (let dx = -r; dx <= r; dx++) {",
          "        const testX = x + dx;",
          "        const testY = y + dy;",
          "        if (testX >= 0 && testX < GRID_WIDTH && testY >= 0 && testY < GRID_HEIGHT) {",
          "          const testIndex = testY * GRID_WIDTH + testX;",
          "          if (!isSolid[testIndex]) return testIndex;",
          "        }",
          "      }",
          "    }",
          "  }",
          "  return 0;",
          "}",
        ],
        explanations: [
          "Integrate forces around cylinder surface",
          "Calculate drag and lift coefficients",
          "Compute Reynolds number for flow regime",
          "Compare with experimental correlations",
          "Enable engineering design decisions",
        ],
      },
      executionResult: {
        description:
          "Complete engineering analysis - drag, lift, and flow regime identified!",
        visualEffect: "analysis",
      },
    },
  ];

  // Simple simulation for visualization with proper flow physics
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d")!;
    canvas.width = 800;
    canvas.height = 400;

    // Simulation parameters
    const GRID_W = 100;
    const GRID_H = 50;
    const CELL_SIZE = canvas.width / GRID_W;

    // Flow field arrays
    let u = new Float32Array(GRID_W * GRID_H);
    let v = new Float32Array(GRID_W * GRID_H);
    let u_prev = new Float32Array(GRID_W * GRID_H);
    let v_prev = new Float32Array(GRID_W * GRID_H);
    let obstacle = new Array(GRID_W * GRID_H).fill(false);

    // Obstacle parameters
    const obX = Math.floor(GRID_W * 0.3);
    const obY = Math.floor(GRID_H * 0.5);
    const obR = 6;

    // Initialize flow field and obstacle
    const initializeSimulation = () => {
      // Set initial uniform flow
      for (let i = 0; i < GRID_W * GRID_H; i++) {
        u[i] = currentStep >= 1 ? 3.0 : 0.0;
        v[i] = 0.0;
      }

      // Create obstacle
      if (currentStep >= 2) {
        for (let y = 0; y < GRID_H; y++) {
          for (let x = 0; x < GRID_W; x++) {
            const dx = x - obX;
            const dy = y - obY;
            const dist = Math.sqrt(dx * dx + dy * dy);
            const idx = y * GRID_W + x;
            if (dist < obR) {
              obstacle[idx] = true;
              u[idx] = 0;
              v[idx] = 0;
            }
          }
        }
      }
    };

    // Boundary conditions
    const setBoundaryConditions = () => {
      // Inlet (left boundary)
      if (currentStep >= 1) {
        for (let y = 0; y < GRID_H; y++) {
          const idx = y * GRID_W + 0;
          u[idx] = 3.0;
          v[idx] = 0.0;
        }
      }

      // Walls (top and bottom)
      for (let x = 0; x < GRID_W; x++) {
        u[x] = 0; // top wall
        v[x] = 0;
        const bottomIdx = (GRID_H - 1) * GRID_W + x;
        u[bottomIdx] = 0; // bottom wall
        v[bottomIdx] = 0;
      }

      // Obstacle boundaries
      if (currentStep >= 2) {
        for (let i = 0; i < GRID_W * GRID_H; i++) {
          if (obstacle[i]) {
            u[i] = 0;
            v[i] = 0;
          }
        }
      }
    };

    // Simple diffusion step
    const diffuse = (
      field: Float32Array,
      prev: Float32Array,
      viscosity: number,
      dt: number
    ) => {
      if (currentStep < 6) return;

      const a = dt * viscosity * GRID_W * GRID_H;
      const c = 1 + 4 * a;

      for (let iter = 0; iter < 10; iter++) {
        for (let y = 1; y < GRID_H - 1; y++) {
          for (let x = 1; x < GRID_W - 1; x++) {
            const idx = y * GRID_W + x;
            if (obstacle[idx]) continue;

            field[idx] =
              (prev[idx] +
                a *
                  (field[idx - 1] +
                    field[idx + 1] +
                    field[idx - GRID_W] +
                    field[idx + GRID_W])) /
              c;
          }
        }
      }
    };

    // Advection step
    const advect = (
      field: Float32Array,
      prev: Float32Array,
      u_vel: Float32Array,
      v_vel: Float32Array,
      dt: number
    ) => {
      if (currentStep < 4) return;

      for (let y = 1; y < GRID_H - 1; y++) {
        for (let x = 1; x < GRID_W - 1; x++) {
          const idx = y * GRID_W + x;
          if (obstacle[idx]) continue;

          const backX = x - dt * u_vel[idx];
          const backY = y - dt * v_vel[idx];

          const x0 = Math.max(0, Math.min(GRID_W - 1, Math.floor(backX)));
          const y0 = Math.max(0, Math.min(GRID_H - 1, Math.floor(backY)));
          const x1 = Math.min(GRID_W - 1, x0 + 1);
          const y1 = Math.min(GRID_H - 1, y0 + 1);

          const fx = backX - x0;
          const fy = backY - y0;

          field[idx] =
            prev[y0 * GRID_W + x0] * (1 - fx) * (1 - fy) +
            prev[y0 * GRID_W + x1] * fx * (1 - fy) +
            prev[y1 * GRID_W + x0] * (1 - fx) * fy +
            prev[y1 * GRID_W + x1] * fx * fy;
        }
      }
    };

    // Pressure projection
    const project = () => {
      if (currentStep < 5) return;

      let divergence = new Float32Array(GRID_W * GRID_H);
      let pressure = new Float32Array(GRID_W * GRID_H);

      // Calculate divergence
      for (let y = 1; y < GRID_H - 1; y++) {
        for (let x = 1; x < GRID_W - 1; x++) {
          const idx = y * GRID_W + x;
          if (obstacle[idx]) continue;

          divergence[idx] =
            0.5 * (u[idx + 1] - u[idx - 1] + v[idx + GRID_W] - v[idx - GRID_W]);
        }
      }

      // Solve pressure Poisson equation
      for (let iter = 0; iter < 20; iter++) {
        for (let y = 1; y < GRID_H - 1; y++) {
          for (let x = 1; x < GRID_W - 1; x++) {
            const idx = y * GRID_W + x;
            if (obstacle[idx]) continue;

            pressure[idx] =
              0.25 *
              (pressure[idx - 1] +
                pressure[idx + 1] +
                pressure[idx - GRID_W] +
                pressure[idx + GRID_W] -
                divergence[idx]);
          }
        }
      }

      // Apply pressure correction
      for (let y = 1; y < GRID_H - 1; y++) {
        for (let x = 1; x < GRID_W - 1; x++) {
          const idx = y * GRID_W + x;
          if (obstacle[idx]) continue;

          u[idx] -= 0.5 * (pressure[idx + 1] - pressure[idx - 1]);
          v[idx] -= 0.5 * (pressure[idx + GRID_W] - pressure[idx - GRID_W]);
        }
      }
    };

    // Simulation step
    const simulationStep = () => {
      if (!isRunning || currentStep < 1) return;

      const dt = 0.1;
      const viscosity = 0.0001;

      // Copy current to previous
      u_prev.set(u);
      v_prev.set(v);

      // Advection
      if (currentStep >= 4) {
        advect(u, u_prev, u_prev, v_prev, dt);
        advect(v, v_prev, u_prev, v_prev, dt);
        setBoundaryConditions();
      }

      // Diffusion
      if (currentStep >= 6) {
        u_prev.set(u);
        v_prev.set(v);
        diffuse(u, u_prev, viscosity, dt);
        diffuse(v, v_prev, viscosity, dt);
        setBoundaryConditions();
      }

      // Projection
      if (currentStep >= 5) {
        project();
        setBoundaryConditions();
      }
    };

    const render = () => {
      // Clear canvas
      ctx.fillStyle = "#1f2937";
      ctx.fillRect(0, 0, canvas.width, canvas.height);

      // Draw grid (step 0+)
      if (currentStep >= 0) {
        ctx.strokeStyle = "#374151";
        ctx.lineWidth = 0.5;
        for (let x = 0; x < canvas.width; x += CELL_SIZE) {
          ctx.beginPath();
          ctx.moveTo(x, 0);
          ctx.lineTo(x, canvas.height);
          ctx.stroke();
        }
        for (let y = 0; y < canvas.height; y += CELL_SIZE) {
          ctx.beginPath();
          ctx.moveTo(0, y);
          ctx.lineTo(canvas.width, y);
          ctx.stroke();
        }

        // Add progression animation - fill grid from left to right
        const progressWidth = (canvas.width * (currentStep + 1)) / steps.length;
        const gradient = ctx.createLinearGradient(0, 0, progressWidth, 0);
        gradient.addColorStop(0, "rgba(59, 130, 246, 0.1)");
        gradient.addColorStop(1, "rgba(59, 130, 246, 0.05)");
        ctx.fillStyle = gradient;
        ctx.fillRect(0, 0, progressWidth, canvas.height);
      }

      // Draw velocity field as arrows (step 1+)
      if (currentStep >= 1) {
        ctx.strokeStyle = "#60a5fa";
        ctx.fillStyle = "#60a5fa";
        ctx.lineWidth = 1;

        const skipCells = Math.max(1, Math.floor(8 / CELL_SIZE));
        for (let y = skipCells; y < GRID_H; y += skipCells) {
          for (let x = skipCells; x < GRID_W; x += skipCells) {
            const idx = y * GRID_W + x;
            if (obstacle[idx]) continue;

            const centerX = x * CELL_SIZE + CELL_SIZE / 2;
            const centerY = y * CELL_SIZE + CELL_SIZE / 2;

            const vel_magnitude = Math.sqrt(u[idx] * u[idx] + v[idx] * v[idx]);
            if (vel_magnitude < 0.1) continue;

            const arrow_length = Math.min(CELL_SIZE * 2, vel_magnitude * 8);
            const angle = Math.atan2(v[idx], u[idx]);

            const endX = centerX + arrow_length * Math.cos(angle);
            const endY = centerY + arrow_length * Math.sin(angle);

            // Draw arrow shaft
            ctx.beginPath();
            ctx.moveTo(centerX, centerY);
            ctx.lineTo(endX, endY);
            ctx.stroke();

            // Draw arrow head
            const headSize = 4;
            const headAngle1 = angle + Math.PI * 0.8;
            const headAngle2 = angle - Math.PI * 0.8;

            ctx.beginPath();
            ctx.moveTo(endX, endY);
            ctx.lineTo(
              endX + headSize * Math.cos(headAngle1),
              endY + headSize * Math.sin(headAngle1)
            );
            ctx.lineTo(
              endX + headSize * Math.cos(headAngle2),
              endY + headSize * Math.sin(headAngle2)
            );
            ctx.closePath();
            ctx.fill();
          }
        }
      }

      // Draw obstacle (step 2+)
      if (currentStep >= 2) {
        const centerX = obX * CELL_SIZE + CELL_SIZE / 2;
        const centerY = obY * CELL_SIZE + CELL_SIZE / 2;
        const radius = obR * CELL_SIZE;

        // Draw obstacle with gradient
        const gradient = ctx.createRadialGradient(
          centerX,
          centerY,
          0,
          centerX,
          centerY,
          radius
        );
        gradient.addColorStop(0, "#ef4444");
        gradient.addColorStop(1, "#dc2626");
        ctx.fillStyle = gradient;

        ctx.beginPath();
        ctx.arc(centerX, centerY, radius, 0, 2 * Math.PI);
        ctx.fill();

        // Add border
        ctx.strokeStyle = "#991b1b";
        ctx.lineWidth = 2;
        ctx.stroke();

        // Add label
        ctx.fillStyle = "#ffffff";
        ctx.font = "bold 12px monospace";
        ctx.textAlign = "center";
        ctx.fillText("CYLINDER", centerX, centerY + 4);
      }

      // Draw streamlines or flow visualization (step 3+)
      if (currentStep >= 3 && isRunning) {
        // Draw some streamlines
        ctx.strokeStyle = "rgba(251, 191, 36, 0.8)";
        ctx.lineWidth = 2;

        for (let startY = 10; startY < GRID_H - 10; startY += 8) {
          let x = 5;
          let y = startY;

          ctx.beginPath();
          ctx.moveTo(x * CELL_SIZE, y * CELL_SIZE);

          for (let step = 0; step < 200 && x < GRID_W - 5; step++) {
            const idx = Math.floor(y) * GRID_W + Math.floor(x);
            if (idx < 0 || idx >= GRID_W * GRID_H || obstacle[idx]) break;

            const dt = 0.1;
            x += u[idx] * dt;
            y += v[idx] * dt;

            if (x < 0 || x >= GRID_W || y < 0 || y >= GRID_H) break;

            ctx.lineTo(x * CELL_SIZE, y * CELL_SIZE);
          }

          ctx.stroke();
        }
      }

      // Show pressure field visualization (step 5+)
      if (currentStep >= 5 && showPhysics && isRunning) {
        for (let y = 0; y < GRID_H; y++) {
          for (let x = 0; x < GRID_W; x++) {
            const idx = y * GRID_W + x;
            if (obstacle[idx]) continue;

            const vel_mag = Math.sqrt(u[idx] * u[idx] + v[idx] * v[idx]);
            const pressure_estimate = 1 - 0.5 * vel_mag * vel_mag; // Bernoulli approximation

            const intensity = Math.max(
              0,
              Math.min(1, (pressure_estimate + 1) * 0.5)
            );
            ctx.fillStyle = `rgba(147, 51, 234, ${intensity * 0.3})`;
            ctx.fillRect(x * CELL_SIZE, y * CELL_SIZE, CELL_SIZE, CELL_SIZE);
          }
        }
      }

      // Show vorticity field (step 8+)
      if (currentStep >= 8 && showPhysics && isRunning) {
        // Calculate and display vorticity
        for (let y = 1; y < GRID_H - 1; y++) {
          for (let x = 1; x < GRID_W - 1; x++) {
            const idx = y * GRID_W + x;
            if (obstacle[idx]) continue;

            // Calculate vorticity: ω = ∂v/∂x - ∂u/∂y
            const dvdx = 0.5 * (v[idx + 1] - v[idx - 1]);
            const dudy = 0.5 * (u[idx + GRID_W] - u[idx - GRID_W]);
            const vorticity = dvdx - dudy;

            const intensity = Math.min(1, Math.abs(vorticity) * 20);
            if (intensity > 0.1) {
              // Blue for positive vorticity, red for negative
              if (vorticity > 0) {
                ctx.fillStyle = `rgba(59, 130, 246, ${intensity * 0.7})`;
              } else {
                ctx.fillStyle = `rgba(239, 68, 68, ${intensity * 0.7})`;
              }
              ctx.fillRect(x * CELL_SIZE, y * CELL_SIZE, CELL_SIZE, CELL_SIZE);
            }
          }
        }
      }

      // Show analysis data overlay (step 9+)
      if (currentStep >= 9 && isRunning) {
        // Display Reynolds number and drag coefficient
        ctx.fillStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillRect(10, 10, 200, 80);

        ctx.fillStyle = "#ffffff";
        ctx.font = "14px monospace";
        ctx.fillText("CFD Analysis:", 15, 30);
        ctx.fillText(`Re ≈ ${Math.floor((obR * 2 * 3.0) / 0.0001)}`, 15, 50);
        ctx.fillText("Cd ≈ 1.2 (theory)", 15, 70);
        ctx.fillText("Vortex shedding active", 15, 90);
      }
    };

    const animate = () => {
      simulationStep();
      render();
      if (isRunning || currentStep === 0) {
        animationRef.current = requestAnimationFrame(animate);
      }
    };

    initializeSimulation();
    setBoundaryConditions();
    animate();

    return () => {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  }, [currentStep, isRunning, showPhysics]);

  const nextStep = () => {
    if (currentStep < steps.length - 1) {
      setCurrentStep(currentStep + 1);
    }
  };

  const prevStep = () => {
    if (currentStep > 0) {
      setCurrentStep(currentStep - 1);
    }
  };

  const resetSimulation = () => {
    setCurrentStep(0);
    setIsRunning(false);
  };

  const toggleSimulation = () => {
    setIsRunning(!isRunning);
  };

  const currentStepData = steps[currentStep];

  return (
    <div className="min-h-screen bg-gray-900 text-white">
      <div className="max-w-7xl mx-auto px-4 py-8">
        {/* Header */}
        <div className="text-center mb-8">
          <h1 className="text-3xl sm:text-4xl font-bold mb-4 bg-gradient-to-r from-blue-400 to-cyan-400 bg-clip-text text-transparent">
            Build Your Own CFD Solver
          </h1>
          <p className="text-gray-400 text-base sm:text-lg">
            Learn CFD by building a fluid dynamics simulator step by step
          </p>
          <div className="mt-4 flex flex-wrap justify-center gap-2 sm:gap-4">
            <span className="bg-blue-600/20 text-blue-400 px-3 py-1 rounded-full text-sm">
              Step {currentStep + 1} of {steps.length}
            </span>
            <span className="bg-green-600/20 text-green-400 px-3 py-1 rounded-full text-sm">
              Interactive Learning
            </span>
            <span className="bg-purple-600/20 text-purple-400 px-3 py-1 rounded-full text-sm">
              {Math.round(((currentStep + 1) / steps.length) * 100)}% Complete
            </span>
          </div>
        </div>{" "}
        {/* Main Layout */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 xl:gap-8">
          {/* Left Column - Visualization */}
          <div className="space-y-6">
            {/* Canvas */}
            <div className="bg-gray-800 rounded-lg p-6">
              <div className="flex flex-wrap gap-3 justify-between items-center mb-4">
                <h2 className="text-xl font-semibold flex items-center gap-2">
                  <Eye className="w-5 h-5" />
                  Live Simulation
                </h2>
                <div className="flex gap-2 flex-wrap">
                  <button
                    onClick={toggleSimulation}
                    className={`px-4 py-2 rounded-lg flex items-center gap-2 transition-colors ${
                      isRunning
                        ? "bg-red-600 hover:bg-red-700"
                        : "bg-green-600 hover:bg-green-700"
                    }`}
                  >
                    {isRunning ? (
                      <Play className="w-4 h-4" />
                    ) : (
                      <Play className="w-4 h-4" />
                    )}
                    {isRunning ? "Running" : "Start"}
                  </button>
                  <button
                    onClick={resetSimulation}
                    className="px-4 py-2 bg-gray-600 hover:bg-gray-700 rounded-lg flex items-center gap-2"
                  >
                    <RotateCcw className="w-4 h-4" />
                    Reset
                  </button>
                </div>
              </div>
              <canvas
                ref={canvasRef}
                className="w-full border border-gray-600 rounded-lg"
                style={{ aspectRatio: "2/1" }}
              />
              <div className="mt-4 text-sm text-gray-400">
                {currentStepData.executionResult.description}
              </div>
            </div>

            {/* Controls */}
            <div className="bg-gray-800 rounded-lg p-6">
              <div className="flex flex-wrap gap-3 justify-between items-center mb-4">
                <button
                  onClick={prevStep}
                  disabled={currentStep === 0}
                  className="px-6 py-3 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-600 disabled:cursor-not-allowed rounded-lg flex items-center gap-2"
                >
                  <ChevronLeft className="w-4 h-4" />
                  Previous
                </button>

                <div className="text-center">
                  <h3 className="text-lg font-semibold">
                    {currentStepData.title}
                  </h3>
                  <p className="text-gray-400 text-sm">
                    {currentStepData.description}
                  </p>
                </div>

                <button
                  onClick={nextStep}
                  disabled={currentStep === steps.length - 1}
                  className="px-6 py-3 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-600 disabled:cursor-not-allowed rounded-lg flex items-center gap-2"
                >
                  Next
                  <ChevronRight className="w-4 h-4" />
                </button>
              </div>

              {/* Progress Bar */}
              <div className="w-full bg-gray-700 rounded-full h-2 mb-4">
                <div
                  className="bg-blue-600 h-2 rounded-full transition-all duration-300"
                  style={{
                    width: `${((currentStep + 1) / steps.length) * 100}%`,
                  }}
                ></div>
              </div>

              {/* Step Indicators */}
              <div className="flex justify-between mb-4">
                {steps.map((step, index) => (
                  <div key={index} className="flex flex-col items-center">
                    <div
                      className={`w-3 h-3 rounded-full transition-colors ${
                        index <= currentStep ? "bg-blue-500" : "bg-gray-600"
                      }`}
                    />
                    <span
                      className={`text-xs mt-1 transition-colors ${
                        index <= currentStep ? "text-blue-400" : "text-gray-500"
                      }`}
                    >
                      {index + 1}
                    </span>
                  </div>
                ))}
              </div>
            </div>
          </div>

          {/* Right Column - Content */}
          <div className="space-y-6">
            {/* Toggle Buttons */}
            <div className="flex gap-2 flex-wrap">
              <button
                onClick={() => setShowPhysics(!showPhysics)}
                className={`px-4 py-2 rounded-lg flex items-center gap-2 ${
                  showPhysics ? "bg-purple-600" : "bg-gray-600"
                }`}
              >
                <Calculator className="w-4 h-4" />
                Physics
              </button>
              <button
                onClick={() => setShowCode(!showCode)}
                className={`px-4 py-2 rounded-lg flex items-center gap-2 ${
                  showCode ? "bg-green-600" : "bg-gray-600"
                }`}
              >
                <Code className="w-4 h-4" />
                Code
              </button>
            </div>

            {/* Physics Explanation */}
            {showPhysics && (
              <div className="bg-purple-900/30 border border-purple-600/30 rounded-lg p-6">
                <h3 className="text-lg font-semibold mb-3 flex items-center gap-2">
                  <BookOpen className="w-5 h-5 text-purple-400" />
                  Physics Behind This Step
                </h3>
                <p className="text-gray-300 leading-relaxed">
                  {currentStepData.physicsExplanation}
                </p>
              </div>
            )}

            {/* Code Changes */}
            {showCode && (
              <div className="bg-gray-800 rounded-lg">
                <div className="px-6 py-4 border-b border-gray-700">
                  <h3 className="text-lg font-semibold flex items-center gap-2">
                    <Code className="w-5 h-5 text-green-400" />
                    Complete CFD Code
                    <span className="text-sm text-gray-400 ml-2">
                      (Current step highlighted)
                    </span>
                  </h3>
                </div>
                <div className="p-4 sm:p-6">
                  <div
                    ref={codeScrollRef}
                    className="bg-gray-900 rounded-lg p-3 sm:p-4 sm:pl-12 max-h-64 sm:max-h-96 overflow-y-auto overflow-x-auto relative"
                  >
                    <pre>
                      <code className="text-xs sm:text-sm">
                        {getCumulativeCode().map((lineObj, index) => (
                          <div
                            key={index}
                            data-is-new={lineObj.isNew}
                            data-step-id={lineObj.stepId}
                            className={`relative ${
                              lineObj.text.trim() === "" ? "h-4" : ""
                            } ${
                              lineObj.isNew
                                ? `bg-green-900/30 border-l-2 border-green-400 pl-3 ${
                                    isScrolling ? "animate-pulse" : ""
                                  }`
                                : ""
                            } ${
                              lineObj.isStepSeparator
                                ? "border-t border-gray-600 pt-2 mt-2"
                                : ""
                            }`}
                          >
                            {/* Line number */}
                            <span className="absolute left-[-35px] text-xs text-gray-600 select-none w-6 text-right">
                              {lineObj.text.trim() !== "" &&
                              !lineObj.isStepSeparator
                                ? getCumulativeCode()
                                    .slice(0, index + 1)
                                    .filter(
                                      (l) =>
                                        l.text.trim() !== "" &&
                                        !l.isStepSeparator
                                    ).length
                                : ""}
                            </span>

                            {lineObj.isStepSeparator ? (
                              <span className="text-cyan-400 font-semibold bg-gray-800/50 px-2 py-1 rounded">
                                {lineObj.text}
                              </span>
                            ) : lineObj.text.startsWith("//") &&
                              !lineObj.isStepSeparator ? (
                              <span className="text-gray-500 italic">
                                <span className="text-gray-600"> </span>
                                {lineObj.text}
                              </span>
                            ) : (
                              <>
                                <span
                                  className={`${
                                    lineObj.isNew
                                      ? "text-green-400 font-semibold"
                                      : "text-gray-600"
                                  }`}
                                >
                                  {lineObj.isNew ? "+" : " "}
                                </span>
                                <span
                                  className={`ml-1 ${
                                    lineObj.isNew
                                      ? "text-gray-100 font-medium"
                                      : "text-gray-500"
                                  }`}
                                >
                                  {lineObj.text}
                                </span>
                              </>
                            )}
                          </div>
                        ))}
                      </code>
                    </pre>

                    {/* Scroll indicator */}
                    <div
                      className={`absolute top-2 right-2 text-xs px-2 py-1 rounded transition-all duration-300 ${
                        isScrolling
                          ? "bg-green-600/80 text-white"
                          : "bg-gray-800/80 text-gray-500"
                      }`}
                    >
                      {isScrolling
                        ? "✨ Scrolling to new code..."
                        : "Scroll to see all code"}
                    </div>
                  </div>
                  <div className="mt-4">
                    <div className="flex justify-between items-start">
                      <div>
                        <h4 className="text-sm font-semibold text-yellow-400 mb-2">
                          Current Step - Key Concepts:
                        </h4>
                        <ul className="text-sm text-gray-400 space-y-1">
                          {currentStepData.codeChanges.explanations.map(
                            (explanation, index) => (
                              <li
                                key={index}
                                className="flex items-start gap-2"
                              >
                                <span className="text-yellow-400 mt-1">•</span>
                                {explanation}
                              </li>
                            )
                          )}
                        </ul>
                      </div>
                      <div className="text-right text-xs text-gray-500">
                        <div>
                          Total lines:{" "}
                          {
                            getCumulativeCode().filter(
                              (l) => l.text.trim() !== ""
                            ).length
                          }
                        </div>
                        <div>
                          New in this step:{" "}
                          {
                            currentStepData.codeChanges.additions.filter(
                              (l) => l.trim() !== ""
                            ).length
                          }
                        </div>
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>
        {/* Footer */}
        <div className="mt-12 text-center">
          <div className="inline-flex items-center gap-2 bg-gray-800 px-4 py-2 rounded-lg">
            <Zap className="w-4 h-4 text-yellow-400" />
            <span className="text-sm text-gray-300">
              Interactive CFD Learning Experience
            </span>
          </div>
        </div>
      </div>
    </div>
  );
};

export default InteractiveCFDBuilder;
